"""
Tests for the Somalier flags report, run entirely offline against the mock fixtures.

The interesting behaviour is all in the grouping: which family a flag lands under, what happens
when a pedigree pair straddles two families, and how a family that has both active and resolved
flags gets split across the two sections. Everything here is a pure function, so no Metamist.
"""

import pytest
from fixtures.somalier_flags import (
    MOCK_ALL_CLEAR_SEQUENCING_GROUPS,
    MOCK_SEQUENCING_GROUPS,
    MOCK_SG_INFOS,
    pedigree_flag,
    sex_flag,
)

from rd_qc.scripts.somalier_flags_report import (
    INLINE_FLAG_LIMIT,
    SgFlags,
    SGInfo,
    _extract_reads,
    _fmt_num,
    _sg_id_rank,
    collect_somalier_flags,
    flag_sg_key,
    group_by_family,
    referenced_sg_ids,
    render_report,
    split_active_resolved,
    split_by_impact,
    summarise_flags,
    summary_message_text,
)
from rd_qc.utils import (
    UNSPECIFIED_RELATED,
    SomalierRelatednessFlag,
    SomalierSelfRelatednessFlag,
    SomalierSexInferenceFlag,
)

# The glance line leads with participants; CPG004 is PID_C in FAM02, CPG010 is PID_H in FAM07.
CROSS_FAMILY_SUBJECT = 'PID_C ↔ PID_H'


def run_pipeline(sequencing_groups=MOCK_SEQUENCING_GROUPS, infos=MOCK_SG_INFOS):
    """The whole report pipeline, as main() runs it but without the two Metamist queries."""
    flagged = [sf for sf in collect_somalier_flags(sequencing_groups) if sf.flags]
    groups = group_by_family(flagged, infos)
    active, resolved = split_active_resolved(groups, infos)
    summary = summarise_flags(flagged, total_sgs=len(sequencing_groups), families_affected=len(active))
    return flagged, active, resolved, summary


def group_by_label(groups) -> dict:
    return {group.label: group for group in groups}


# ---------------------------------------------------------------------------
# Reading flags out of meta
# ---------------------------------------------------------------------------
def test_collect_dispatches_each_category_to_its_dataclass():
    flagged = {sf.sg_id: sf for sf in collect_somalier_flags(MOCK_SEQUENCING_GROUPS)}

    assert isinstance(flagged['CPG001'].flags[0], SomalierSexInferenceFlag)
    assert isinstance(flagged['CPG002'].flags[0], SomalierSelfRelatednessFlag)
    assert isinstance(flagged['CPG004'].flags[0], SomalierRelatednessFlag)


def test_collect_skips_unrecognised_category_without_raising():
    groups = [{'id': 'CPG001', 'meta': {'somalier_flags': [{'category': 'something_new', 'resolved': False}]}}]

    collected = collect_somalier_flags(groups)

    assert collected[0].flags == ()


def test_collect_skips_malformed_flag_without_raising():
    # A sex flag missing every one of its required measured fields.
    groups = [{'id': 'CPG001', 'meta': {'somalier_flags': [{'category': 'sex_inference_mismatch'}]}}]

    collected = collect_somalier_flags(groups)

    assert collected[0].flags == ()


def test_collect_handles_sg_with_no_meta_at_all():
    collected = collect_somalier_flags([{'id': 'CPG001', 'meta': None}])

    assert collected == [SgFlags(sg_id='CPG001', flags=())]


# ---------------------------------------------------------------------------
# The SG key, and the fallback for flags recorded before it existed
# ---------------------------------------------------------------------------
def test_flag_sg_key_uses_the_recorded_key():
    flag = SomalierRelatednessFlag(**pedigree_flag('CPG004', 'CPG005', 'FAM02', 'parent-child', 'unrelated'))

    assert flag_sg_key(flag, owning_sg_id='CPG004') == 'CPG004_CPG005'


def test_flag_sg_key_falls_back_to_the_pair_ids_for_legacy_flags():
    raw = pedigree_flag('CPG004', 'CPG005', 'FAM02', 'parent-child', 'unrelated')
    flag = SomalierRelatednessFlag(**{**raw, 'sequencing_group_key': ''})

    assert flag_sg_key(flag, owning_sg_id='CPG004') == 'CPG004_CPG005'


def test_flag_sg_key_falls_back_to_the_owning_sg_for_legacy_sex_flags():
    flag = SomalierSexInferenceFlag(**{**sex_flag('CPG001', 'female', 'male'), 'sequencing_group_key': ''})

    assert flag_sg_key(flag, owning_sg_id='CPG001') == 'CPG001'


def test_referenced_sg_ids_includes_the_far_member_of_a_pair():
    # CPG010 and CPG011 only ever appear as the second member of a pair, never as an owning SG.
    flagged, _, _, _ = run_pipeline()

    referenced = referenced_sg_ids(flagged)

    assert 'CPG010' in referenced
    assert 'CPG011' in referenced


# ---------------------------------------------------------------------------
# Grouping
# ---------------------------------------------------------------------------
def test_cross_family_pedigree_flag_lands_in_both_families():
    _, active, _, _ = run_pipeline()
    by_label = group_by_label(active)

    fam02_subjects = [f.subject for f in by_label['FAM02'].flags]
    fam07_subjects = [f.subject for f in by_label['FAM07'].flags]

    assert CROSS_FAMILY_SUBJECT in fam02_subjects
    assert CROSS_FAMILY_SUBJECT in fam07_subjects


def test_cross_family_row_names_the_other_family_in_each_copy():
    _, active, _, _ = run_pipeline()
    by_label = group_by_label(active)

    fam02_row = next(f for f in by_label['FAM02'].flags if f.subject == CROSS_FAMILY_SUBJECT)
    fam07_row = next(f for f in by_label['FAM07'].flags if f.subject == CROSS_FAMILY_SUBJECT)

    assert fam02_row.cross_family == 'FAM07'
    assert fam07_row.cross_family == 'FAM02'


def test_cross_family_flag_is_counted_once_despite_appearing_twice():
    _, active, _, summary = run_pipeline()

    rendered_rows = sum(len(group.flags) for group in active)

    # Eight rows on the page, but only seven distinct flags: the cross-family pair is shown twice.
    assert rendered_rows == 8
    assert summary['active_flags'] == 7


def test_same_family_pedigree_flag_lands_in_one_group_only():
    _, active, _, _ = run_pipeline()
    same_family_subject = 'PID_C ↔ PID_D'

    holders = [group.label for group in active if any(f.subject == same_family_subject for f in group.flags)]

    assert holders == ['FAM02']


def test_self_relatedness_groups_under_the_family_from_metamist():
    # The flag records participant PID_B but no family; FAM03 comes from the SG info lookup.
    _, active, _, _ = run_pipeline()
    by_label = group_by_label(active)

    self_rows = [f for f in by_label['FAM03'].flags if f.category_key == 'self']

    # Both SGs belong to PID_B, so the subject is the one participant rather than a pair.
    assert [f.subject for f in self_rows] == ['PID_B']


def test_self_relatedness_falls_back_to_a_participant_group_with_no_family():
    _, active, _, _ = run_pipeline()

    participant_groups = [group for group in active if group.key == 'participant:PID_G']

    assert len(participant_groups) == 1
    assert participant_groups[0].label == '(no family) · PID_G'
    assert [f.subject for f in participant_groups[0].flags] == ['PID_G']


def test_group_carries_the_sg_info_for_both_members_of_a_pair():
    _, active, _, _ = run_pipeline()
    by_label = group_by_label(active)

    sg_ids = [info.sg_id for info in by_label['FAM02'].sg_infos]

    assert set(sg_ids) >= {'CPG004', 'CPG005', 'CPG010'}


def test_groups_sort_worst_first():
    _, active, _, _ = run_pipeline()

    totals = [sum(group.counts.values()) for group in active]

    assert totals == sorted(totals, reverse=True)


# ---------------------------------------------------------------------------
# Active / resolved split
# ---------------------------------------------------------------------------
def test_family_with_both_kinds_appears_in_both_sections():
    _, active, resolved, _ = run_pipeline()

    assert 'FAM02' in group_by_label(active)
    assert 'FAM02' in group_by_label(resolved)


def test_each_section_shows_only_its_own_flags():
    _, active, resolved, _ = run_pipeline()

    assert all(not f.resolved for group in active for f in group.flags)
    assert all(f.resolved for group in resolved for f in group.flags)


def test_counts_are_recomputed_per_section():
    _, active, resolved, _ = run_pipeline()

    # FAM02 has two active pedigree flags, and one resolved pedigree plus one resolved sex.
    assert group_by_label(active)['FAM02'].counts == {'pedigree': 2}
    assert group_by_label(resolved)['FAM02'].counts == {'sex': 1, 'pedigree': 1}


def test_summary_counts():
    _, active, _, summary = run_pipeline()

    assert summary == {
        'total_sgs': 13,
        'active_flags': 7,
        'active_by_category': {'sex': 2, 'self': 2, 'pedigree': 3},
        # Six of the seven are real disagreements; FAM08's siblings/full-siblings pair is not.
        # That one is a pedigree flag, so pedigree drops to 2 once refinements come out.
        'active_conflicts_by_category': {'sex': 2, 'self': 2, 'pedigree': 2},
        'active_conflicts': 6,
        'active_refinements': 1,
        'families_affected': len(active),
        'resolved_flags': 2,
    }


# ---------------------------------------------------------------------------
# Number formatting
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    ('value', 'expected'),
    [
        (0.616722, '0.62'),
        (0.0616722, '0.062'),
        (64.579124, '64.58'),
        (4821.0, '4821'),
        (4821, '4821'),
        (1.5, '1.5'),
        (0.9, '0.9'),
        (2.0, '2'),
        (None, '—'),
        ('unknown', 'unknown'),
    ],
)
def test_fmt_num(value, expected):
    assert _fmt_num(value) == expected


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------
def render_fixture_html(sequencing_groups=MOCK_SEQUENCING_GROUPS) -> str:
    _, active, resolved, summary = run_pipeline(sequencing_groups)
    return render_report('mock-dataset', active, resolved, summary=summary, generated_at='2026-09-14T10:00:00+00:00')


def test_render_leaves_no_unresolved_template_syntax():
    html = render_fixture_html()

    assert '{{' not in html
    assert '{%' not in html
    assert 'Undefined' not in html


@pytest.mark.parametrize(
    'expected',
    [
        'FAM02',
        'FAM07',
        '(no family)',
        'PID_G',
        'cross-family',
        'expected parent-child / measured unrelated',
        'provided unknown / inferred male',
        'relatedness 0.62 (expected ~1.0, threshold 0.9)',
        'Resolved &mdash; past incidents',
        'Families affected',
    ],
)
def test_render_includes_expected_content(expected):
    assert expected in render_fixture_html()


@pytest.mark.parametrize('removed', ['badge-fail', 'badge-warn', 'active_cram', 'badge-source'])
def test_render_has_no_severity_or_source_leftovers(removed):
    assert removed not in render_fixture_html()


def test_all_clear_banner_when_nothing_is_active():
    html = render_fixture_html(MOCK_ALL_CLEAR_SEQUENCING_GROUPS)

    assert 'All clear' in html
    assert 'Pedigree conflicts' in html


def test_resolved_section_is_absent_when_there_is_nothing_resolved():
    # Both the 'past incidents' card subtext and the .resolved-section CSS rule are always present,
    # so this has to look for the section's own heading and wrapper markup.
    html = render_fixture_html(MOCK_ALL_CLEAR_SEQUENCING_GROUPS)

    assert 'Resolved &mdash; past incidents' not in html
    assert 'class="resolved-section"' not in html


def test_inline_flag_lines_are_capped_with_a_more_link():
    # A real family can carry 50+ pedigree mismatches; rendering them all inline buries the rest.
    many = [
        pedigree_flag(f'CPG{i:03d}', f'CPG{i + 1:03d}', 'FAM09', expected='unrelated', inferred='full siblings')
        for i in range(1, INLINE_FLAG_LIMIT + 4)
    ]
    groups = [{'id': 'CPG001', 'meta': {'somalier_flags': many}}]
    infos = {
        f'CPG{i:03d}': SGInfo(
            **{**vars(MOCK_SG_INFOS['CPG004']), 'sg_id': f'CPG{i:03d}', 'family_external_id': 'FAM09'}
        )
        for i in range(1, INLINE_FLAG_LIMIT + 5)
    }
    flagged = [sf for sf in collect_somalier_flags(groups) if sf.flags]
    active, resolved = split_active_resolved(group_by_family(flagged, infos), infos)
    summary = summarise_flags(flagged, total_sgs=1, families_affected=len(active))
    html = render_report('mock', active, resolved, summary=summary, generated_at='2026-09-14T10:00:00+00:00')

    assert active[0].total == INLINE_FLAG_LIMIT + 3
    assert '+3 more &mdash; click to expand' in html
    # Every flag still reaches the page, just via the expanded detail table.
    assert html.count('full siblings') > INLINE_FLAG_LIMIT


def test_split_by_impact_separates_conflicts_from_refinements():
    _, active, _, _ = run_pipeline()

    conflicts, refinements = split_by_impact(active, {info.sg_id: info for g in active for info in g.sg_infos})

    # FAM08's only flag is siblings -> full siblings, so it appears solely in the refinements side.
    assert 'FAM08' not in group_by_label(conflicts)
    assert 'FAM08' in group_by_label(refinements)
    assert all(f.impact == 'conflict' for g in conflicts for f in g.flags)
    assert all(f.impact == 'refinement' for g in refinements for f in g.flags)


def test_refinements_render_in_their_own_section_with_the_explanation():
    html = render_fixture_html()

    assert 'Pedigree refinements' in html
    assert 'Pedigree conflicts' in html
    # The section explains the missing-parent cause, which is why these rows exist at all.
    # Matched within one line, since the template's prose is hard-wrapped.
    assert 'whenever a parent is missing' in html


def test_a_dataset_of_only_refinements_still_shows_the_all_clear():
    # No conflicts means nothing needs a decision, even though flags exist.
    only_refinement = [{'id': 'CPG012', 'meta': MOCK_SEQUENCING_GROUPS[11]['meta']}]
    html = render_fixture_html(only_refinement)

    assert 'All clear' in html
    assert 'Pedigree refinements' in html


# ---------------------------------------------------------------------------
# Slack summary message
# ---------------------------------------------------------------------------
FLAGS_URL = 'https://main-web.populationgenomics.org.au/mock/flags.html'
SOMALIER_URL = 'https://main-web.populationgenomics.org.au/mock/somalier.html'


def message_for(summary, previous_analysis=None) -> str:
    return summary_message_text(
        'mock-dataset',
        flags_html_url=FLAGS_URL,
        somalier_html_url=SOMALIER_URL,
        seq_type='genome',
        seq_tech='short-read',
        summary=summary,
        previous_analysis=previous_analysis,
    )


def previous_report(summary: dict, completed: str = '2026-08-20T03:00:00+00:00') -> dict:
    return {'id': 99, 'timestampCompleted': completed, 'meta': {'summary': summary}}


def test_message_leads_with_the_dataset_and_both_report_links():
    _, _, _, summary = run_pipeline()

    lines = message_for(summary).splitlines()

    assert lines[0] == f'*[mock-dataset]* <{FLAGS_URL}|Somalier flags report (genome | short-read)>'
    assert SOMALIER_URL in lines[1]


def test_message_counts_conflicts_not_total_flags():
    _, active, _, summary = run_pipeline()

    text = message_for(summary)

    # 7 active flags, but only 6 are conflicts. The headline must not claim 7.
    assert '*6 active conflicts*' in text
    assert f'{len(active)} families' in text
    assert '7 active conflicts' not in text


def test_message_breaks_conflicts_down_by_category():
    _, _, _, summary = run_pipeline()

    text = message_for(summary)

    # The one refinement is a pedigree flag, so pedigree reads 2 here and not 3.
    assert ' - Sex inference: 2' in text
    assert ' - Self-relatedness: 2' in text
    assert ' - Pedigree relatedness: 2' in text


def test_message_reports_refinements_separately_from_conflicts():
    _, _, _, summary = run_pipeline()

    text = message_for(summary)

    assert '1 pedigree refinement' in text
    assert 'less specific' in text


def test_message_is_an_all_clear_when_nothing_is_active():
    _, _, _, summary = run_pipeline(MOCK_ALL_CLEAR_SEQUENCING_GROUPS)

    text = message_for(summary)

    assert '✅' in text
    assert 'No active Somalier flags' in text
    assert 'conflicts*' not in text


def test_message_is_an_all_clear_when_only_refinements_are_active():
    summary = {
        'total_sgs': 20,
        'active_flags': 3,
        'active_by_category': {'sex': 0, 'self': 0, 'pedigree': 3},
        'active_conflicts_by_category': {'sex': 0, 'self': 0, 'pedigree': 0},
        'active_conflicts': 0,
        'active_refinements': 3,
        'families_affected': 2,
        'resolved_flags': 0,
    }

    text = message_for(summary)

    # Nothing needs a decision, but the refinements still get their line.
    assert 'No conflicts' in text
    assert '3 pedigree refinements' in text


def test_message_says_nothing_about_changes_without_a_previous_report():
    _, _, _, summary = run_pipeline()

    text = message_for(summary)

    assert 'since the last report' not in text


def test_message_reports_what_is_new_since_the_previous_report():
    _, _, _, summary = run_pipeline()
    previous = previous_report(
        {
            'total_sgs': 11,
            'active_conflicts': 4,
            'active_refinements': 0,
            'families_affected': summary['families_affected'] - 1,
            'resolved_flags': 0,
        }
    )

    text = message_for(summary, previous)

    assert 'Changes since the last report on 2026-08-20' in text
    assert '+2 new sequencing groups' in text
    assert '+1 additional family flagged' in text
    assert '+2 new conflicts' in text
    assert '+1 new pedigree refinement' in text


def test_message_reports_no_change_when_the_summary_is_identical():
    _, _, _, summary = run_pipeline()

    text = message_for(summary, previous_report(dict(summary)))

    assert 'No change since the last report on 2026-08-20' in text
    assert 'Changes since' not in text


def test_message_reports_flags_that_have_been_fixed():
    _, _, _, summary = run_pipeline()
    previous = previous_report({**summary, 'active_conflicts': 9, 'resolved_flags': 0})

    text = message_for(summary, previous)

    assert '3 fewer conflicts' in text
    assert f'{summary["resolved_flags"]} more flags resolved' in text


def test_message_survives_a_previous_summary_from_before_the_conflict_split():
    _, _, _, summary = run_pipeline()
    # Reports registered before conflicts/refinements existed only carry these keys.
    previous = previous_report({'total_sgs': 13, 'active_flags': 7, 'resolved_flags': 2})

    text = message_for(summary, previous)

    assert 'since the last report on 2026-08-20' in text


# ---------------------------------------------------------------------------
# Wording of an unspecified expectation
# ---------------------------------------------------------------------------
def refinement_row():
    """FAM08's only flag: no recorded path between the pair, measured as siblings."""
    _, active, _, _ = run_pipeline()
    return group_by_label(active)['FAM08'].flags[0]


def test_an_unspecified_expectation_reads_as_no_relationship_provided():
    row = refinement_row()

    # 'expected related at unknown level / measured siblings' is accurate but unclear: the point is
    # that the pedigree has no data here, not that it asserted something vague.
    assert row.result == 'No relationship provided / measured siblings'
    assert 'related at unknown level' not in row.result


def test_the_detail_table_uses_the_same_wording():
    assert refinement_row().expected == 'No relationship provided'


def test_a_stated_expectation_is_left_alone():
    _, active, _, _ = run_pipeline()
    row = group_by_label(active)['FAM02'].flags[0]

    assert row.result == 'expected parent-child / measured unrelated'
    assert row.expected == 'parent-child'


def test_relabelling_does_not_touch_the_stored_relationship():
    # The stored value is peddy's vocabulary, pinned in PEDDY_RELATIONSHIPS and part of the flag's
    # reconciliation identity, so only the display may change.
    flagged = {sf.sg_id: sf for sf in collect_somalier_flags(MOCK_SEQUENCING_GROUPS)}

    assert flagged['CPG012'].flags[0].expected_relationship == UNSPECIFIED_RELATED


def test_the_raw_relationship_is_still_searchable():
    # A collaborator pasting the peddy string, or a saved filter, must still find the row.
    assert 'related at unknown level' in refinement_row().search_blob


def test_the_refinements_blurb_describes_the_case_that_actually_occurs():
    html = render_fixture_html()

    # 'siblings' -> 'full siblings' is satisfied by EXPECTED_DEGREES now, so it never reaches the
    # refinements section and must not be described as the common case.
    assert 'full siblings' not in html
    assert 'no relationship' in html.lower()


# ---------------------------------------------------------------------------
# Glance line leads with the identifiers collaborators use
# ---------------------------------------------------------------------------
def test_a_pedigree_pair_leads_with_participants_and_demotes_the_sg_ids():
    _, active, _, _ = run_pipeline()
    row = next(f for f in group_by_label(active)['FAM02'].flags if f.category_key == 'pedigree')

    assert row.subject == 'PID_C ↔ PID_D'
    assert row.subject_detail == 'CPG004 ↔ CPG005'


def test_a_self_relatedness_pair_leads_with_the_one_participant():
    _, active, _, _ = run_pipeline()
    row = next(f for f in group_by_label(active)['FAM03'].flags if f.category_key == 'self')

    assert row.subject == 'PID_B'
    assert row.subject_detail == 'CPG002 ↔ CPG003'


def test_a_pair_with_no_known_participants_falls_back_to_the_sg_ids():
    # Nothing to promote, so the SG ids stay in the lead rather than leaving the line blank.
    flag = SomalierRelatednessFlag(**pedigree_flag('CPG404', 'CPG405', 'FAM99', 'parent-child', 'unrelated'))
    rows = group_by_family([SgFlags(sg_id='CPG404', flags=(flag,))], {})

    assert rows[0].flags[0].subject == 'CPG404 ↔ CPG405'
    assert rows[0].flags[0].subject_detail == ''


# ---------------------------------------------------------------------------
# Read files
# ---------------------------------------------------------------------------
def fastq_assays(reads):
    return [{'meta': {'reads_type': 'fastq', 'reads': reads}}]


def test_fastq_reads_keep_their_pairing_and_carry_size_and_date():
    _, pairs, _ = _extract_reads(
        fastq_assays(
            [
                {
                    'basename': 'S1_R1.fastq.gz',
                    'size': 22280401447,
                    'datetime_added': '2024-06-14T03:27:16.749000+00:00',
                },
                {
                    'basename': 'S1_R2.fastq.gz',
                    'size': 22512341519,
                    'datetime_added': '2024-06-14T03:23:27.163000+00:00',
                },
            ]
        )
    )

    assert len(pairs) == 1
    r1, r2 = pairs[0]
    assert (r1.name, r1.size, r1.date) == ('S1_R1.fastq.gz', '20.75 GiB', '2024-06-14')
    assert (r2.name, r2.size, r2.date) == ('S1_R2.fastq.gz', '20.97 GiB', '2024-06-14')


def test_reads_recorded_without_a_size_or_date_still_appear():
    # Some uploads carry datetime_added: null, so neither field can be assumed present.
    _, pairs, _ = _extract_reads(fastq_assays([{'basename': 'x_R1.fq.gz', 'size': None, 'datetime_added': None}]))

    assert (pairs[0][0].name, pairs[0][0].size, pairs[0][0].date) == ('x_R1.fq.gz', '', '')


def test_an_odd_number_of_fastqs_does_not_invent_a_partner():
    _, pairs, _ = _extract_reads(
        fastq_assays([{'basename': f'S_{n}.fq.gz'} for n in ('1_R1', '1_R2', '2_R1')]),
    )

    assert [[r.name for r in group] for group in pairs] == [['S_1_R1.fq.gz', 'S_1_R2.fq.gz'], ['S_2_R1.fq.gz']]


def test_a_plain_path_string_still_yields_a_name():
    crams, _, _ = _extract_reads([{'meta': {'reads_type': 'cram', 'reads': ['gs://bucket/path/S1.cram']}}])

    assert [c.name for c in crams] == ['S1.cram']


def test_each_read_file_renders_on_its_own_line():
    html = render_fixture_html()

    # Previously R1 and R2 were joined with a slash onto one line, which was unreadable once the
    # real filenames ran to 90 characters.
    assert '<div class="read">' in html
    assert 'EXT_A_R1.fastq.gz</span>' in html
    assert '&nbsp;/&nbsp;' not in html


def test_the_rendered_read_lines_show_size_and_date_when_known():
    html = render_fixture_html()

    # Unlabelled: a GiB figure and a date read as themselves.
    assert '<span class="read-meta">1.50 GiB</span>' in html
    assert '<span class="read-meta">2026-06-01</span>' in html


def test_a_read_with_no_size_or_date_renders_the_name_alone():
    # The fixture's R2 has neither, so it must not emit an empty meta span.
    html = render_fixture_html()
    r2_line = next(line for line in html.splitlines() if 'EXT_A_R2.fastq.gz' in line)

    assert 'read-meta' not in r2_line


# ---------------------------------------------------------------------------
# Sorting by newest sequencing group
# ---------------------------------------------------------------------------
def test_sg_id_rank_is_numeric_not_lexicographic():
    # Real datasets carry both 5- and 6-digit IDs, so a string sort would put CPG99999 first.
    assert _sg_id_rank('CPG100000') > _sg_id_rank('CPG99999')


def test_sg_id_rank_of_an_unparseable_id_sorts_last():
    # Descending sort puts 0 at the bottom, which is where a malformed ID belongs.
    assert _sg_id_rank('not-an-sg') == 0


def test_group_carries_the_rank_of_its_newest_sg():
    _, active, _, _ = run_pipeline()

    # FAM02's flags span CPG004, CPG005 and CPG010, so the newest is CPG010.
    assert group_by_label(active)['FAM02'].newest_sg_rank == 10
