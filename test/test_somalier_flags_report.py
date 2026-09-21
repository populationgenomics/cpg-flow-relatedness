"""
Tests for the Somalier flags report, run entirely offline against the mock fixtures.

The interesting behaviour is all in the grouping: which family a flag lands under, what happens
when a pedigree pair straddles two families, and how a family that has both active and resolved
flags gets split across the two sections. Everything here is a pure function, so no Metamist.
"""

import re

import pytest
from fixtures.somalier_flags import (
    MOCK_ALL_CLEAR_SEQUENCING_GROUPS,
    MOCK_SEQUENCING_GROUPS,
    MOCK_SG_INFOS,
    pedigree_flag,
    sex_flag,
)

from rd_qc.scripts.somalier_flags_report import (
    IMPACT_CONFLICT,
    IMPACT_REFINEMENT,
    IMPACT_SAME_INDIVIDUAL,
    INLINE_FLAG_LIMIT,
    MIN_GROUPS_FOR_FILTER_BAR,
    SgFlags,
    SGInfo,
    _extract_reads,
    _fmt_num,
    _impact_of,
    _is_same_individual,
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
    FIRST_DEGREE_MIN_RELATEDNESS,
    IDENTICAL_MIN_RELATEDNESS,
    RELATEDNESS_BANDS,
    SECOND_DEGREE_MIN_RELATEDNESS,
    THIRD_DEGREE_MIN_RELATEDNESS,
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
    summary = summarise_flags(flagged, total_sgs=len(sequencing_groups), families_affected=len(active), infos=infos)
    return flagged, active, resolved, summary


def group_by_label(groups) -> dict:
    return {group.label: group for group in groups}


def blurb_text(html: str, heading: str) -> str:
    """
    The plain text of the blurb paragraph under a section heading, with markup stripped.

    Scoping to one blurb keeps a section's wording from being satisfied by prose elsewhere on the
    page, and stripping tags keeps these assertions from breaking when a phrase gains or loses an
    <em>. The wording itself is free to change; what these tests pin is that a corrected claim
    does not quietly revert.

    Bounded at the section's own table wrapper rather than running to the next blurb on the page,
    so a section whose blurb was deleted reads as empty instead of silently borrowing the next
    section's paragraph and passing for the wrong reason.
    """
    section = html.split(heading)[1].split('<div data-section=', maxsplit=1)[0]
    if '<p class="blurb">' not in section:
        return ''
    paragraph = section.split('<p class="blurb">')[1].split('</p>')[0]
    return ' '.join(re.sub(r'<[^>]+>', ' ', paragraph).split())


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

    # Nine rows across the active groups, but only eight distinct flags: the cross-family pair is
    # shown twice. FAM09's same-individual pair is one of the nine, since this counts active
    # groups before the impact split.
    assert rendered_rows == 9
    assert summary['active_flags'] == 8


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
        'total_sgs': 15,
        'active_flags': 8,
        'active_by_category': {'sex': 2, 'self': 2, 'pedigree': 4},
        # Six of the eight are real disagreements. FAM08's 'no relationship provided' pair is a
        # refinement, and FAM09's pair is two SGs of one person, so pedigree drops to 2 once both
        # come out.
        'active_conflicts_by_category': {'sex': 2, 'self': 2, 'pedigree': 2},
        'active_conflicts': 6,
        'active_refinements': 1,
        'active_same_individual': 1,
        # Counts families with any active flag, refinement- and same-individual-only included,
        # which is why it tracks len(active) rather than the conflict count.
        'families_affected': len(active),
        'resolved_flags': 2,
    }


def test_every_flag_classifies_into_one_of_the_three_impacts():
    # Fails the moment _impact_of grows a fourth return value, fixture or no fixture, which is
    # what stops active_flags from silently exceeding the sum of its parts.
    flagged, _, _, _ = run_pipeline()

    # The complete lookup, not one rebuilt from the active groups: this classifies every flag
    # including resolved ones, whose SGs need not be referenced by any active flag.
    produced = {_impact_of(flag, MOCK_SG_INFOS) for sf in flagged for flag in sf.flags}

    assert produced <= {IMPACT_CONFLICT, IMPACT_REFINEMENT, IMPACT_SAME_INDIVIDUAL}
    # And the fixture really does exercise all three, so the subset check is not vacuous.
    assert produced == {IMPACT_CONFLICT, IMPACT_REFINEMENT, IMPACT_SAME_INDIVIDUAL}


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
        pedigree_flag(f'CPG{i:03d}', f'CPG{i + 1:03d}', 'FAM_MANY', expected='unrelated', inferred='full siblings')
        for i in range(1, INLINE_FLAG_LIMIT + 4)
    ]
    groups = [{'id': 'CPG001', 'meta': {'somalier_flags': many}}]
    infos = {
        f'CPG{i:03d}': SGInfo(
            **{
                **vars(MOCK_SG_INFOS['CPG004']),
                'sg_id': f'CPG{i:03d}',
                'family_external_id': 'FAM_MANY',
                # Distinct participants: these are meant to be many separate mismatched pairs, not
                # the same-individual case _is_same_individual detects.
                'participant_external_id': f'PID_{i:03d}',
            }
        )
        for i in range(1, INLINE_FLAG_LIMIT + 5)
    }
    flagged = [sf for sf in collect_somalier_flags(groups) if sf.flags]
    active, resolved = split_active_resolved(group_by_family(flagged, infos), infos)
    summary = summarise_flags(flagged, total_sgs=1, families_affected=len(active), infos=infos)
    html = render_report('mock', active, resolved, summary=summary, generated_at='2026-09-14T10:00:00+00:00')

    assert active[0].total == INLINE_FLAG_LIMIT + 3
    assert '+3 more &mdash; click to expand' in html
    # Every flag still reaches the page, just via the expanded detail table.
    assert html.count('full siblings') > INLINE_FLAG_LIMIT


def test_split_by_impact_separates_conflicts_from_refinements():
    _, active, _, _ = run_pipeline()

    conflicts, refinements, _ = split_by_impact(active, {info.sg_id: info for g in active for info in g.sg_infos})

    # FAM08's only flag is 'no relationship provided' -> siblings, so it appears solely in the
    # refinements side.
    assert 'FAM08' not in group_by_label(conflicts)
    assert 'FAM08' in group_by_label(refinements)
    assert all(f.impact == IMPACT_CONFLICT for g in conflicts for f in g.flags)
    assert all(f.impact == IMPACT_REFINEMENT for g in refinements for f in g.flags)


def test_refinements_render_in_their_own_section_with_the_explanation():
    html = render_fixture_html()

    assert 'Pedigree refinements' in html
    assert 'Pedigree conflicts' in html
    # Deliberately not pinning the wording, which is free to change. A section of non-findings
    # with no explanation at all leaves a reader guessing why the rows are there, so only the
    # presence of one is pinned.
    assert blurb_text(html, 'Pedigree refinements') != ''


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

    # 8 active flags, but only 6 are conflicts. The headline must not claim 8.
    assert '*6 active conflicts*' in text
    assert f'{len(active)} families' in text
    assert '8 active conflicts' not in text


def test_message_breaks_conflicts_down_by_category():
    _, _, _, summary = run_pipeline()

    text = message_for(summary)

    # The refinement and the same-individual pair are both pedigree flags, so pedigree reads 2
    # here and not 4.
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
            'total_sgs': 13,
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


def test_the_refinements_blurb_does_not_claim_full_siblings_are_a_refinement():
    # 'siblings' -> 'full siblings' is satisfied by EXPECTED_DEGREES, so that pair never reaches
    # the refinements section and must never be described there as the case that occurs. Scoped
    # to the blurb, because the bands legend names 'full siblings' legitimately.
    blurb = blurb_text(render_fixture_html(), 'Pedigree refinements')

    assert 'full siblings' not in blurb


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


# ---------------------------------------------------------------------------
# Relatedness bands legend
# ---------------------------------------------------------------------------
def test_bands_are_derived_from_the_classifier_thresholds():
    # The legend must quote the numbers infer_degree actually decides on, so a threshold change
    # cannot leave the page describing bands the code no longer uses.
    thresholds = {band.threshold for band in RELATEDNESS_BANDS}

    assert f'>= {IDENTICAL_MIN_RELATEDNESS}' in thresholds
    assert f'>= {FIRST_DEGREE_MIN_RELATEDNESS}' in thresholds
    assert f'>= {SECOND_DEGREE_MIN_RELATEDNESS}' in thresholds
    assert f'>= {THIRD_DEGREE_MIN_RELATEDNESS}' in thresholds
    assert f'< {THIRD_DEGREE_MIN_RELATEDNESS}' in thresholds


def test_every_band_renders_on_the_page():
    html = render_fixture_html()

    assert 'Relatedness bands' in html
    for band in RELATEDNESS_BANDS:
        assert band.label in html
        assert band.expected in html
        assert band.threshold in html
        assert band.examples in html


def test_the_bands_strip_names_the_first_degree_split():
    # 'first-degree' has no DEGREE_* constant because infer_degree splits it on ibs0, so the
    # legend has to spell out both relationships that land in the band.
    html = render_fixture_html()

    assert 'parent-child' in html
    assert 'full siblings' in html


def test_the_refinements_blurb_does_not_promise_every_refinement_is_closable():
    # Many refinements cannot be closed at all: a pedigree encodes relationships only through
    # parent links, so a pair whose connecting individual is absent from the database cannot be
    # stated at all. The blurb must stay hedged rather than promising a pedigree edit resolves
    # them, which is what the original wording did.
    html = render_fixture_html()
    blurb = blurb_text(html, 'Pedigree refinements')

    assert 'the pedigree can be updated to say so' not in html
    assert 'may' in blurb.split()


# ---------------------------------------------------------------------------
# Same individual, multiple sequencing groups
# ---------------------------------------------------------------------------
def test_a_pair_of_sgs_from_one_participant_is_not_a_conflict():
    _, active, _, _ = run_pipeline()
    infos = {info.sg_id: info for g in active for info in g.sg_infos}

    conflicts, refinements, same_individual = split_by_impact(active, infos)

    assert 'FAM09' not in group_by_label(conflicts)
    assert 'FAM09' not in group_by_label(refinements)
    assert 'FAM09' in group_by_label(same_individual)
    # FAM02 carries real conflicts, so an implementation that swept everything into the
    # same-individual bucket would fail here.
    assert 'FAM02' in group_by_label(conflicts)


def test_two_sgs_from_one_participant_are_the_same_individual():
    flag = next(
        f
        for sf in collect_somalier_flags(MOCK_SEQUENCING_GROUPS)
        for f in sf.flags
        if getattr(f, 'sg_id_1', None) == 'CPG014'
    )

    assert _is_same_individual(flag, MOCK_SG_INFOS)


def test_one_resolvable_participant_and_one_missing_is_not_the_same_individual():
    # Asymmetric: CPG014 resolves to PID_M, CPG999 resolves to ''. A set-based rewrite of the
    # guard would behave correctly on the both-empty case and wrongly here.
    flag = SomalierRelatednessFlag(
        category='relatedness_mismatch',
        sg_id_1='CPG014',
        sg_id_2='CPG999',
        family_external_id='FAM09',
        expected_relationship='full siblings',
        inferred_relationship='identical',
        relatedness=0.998,
        ibs0=3,
        ibs2=19871,
    )

    assert not _is_same_individual(flag, MOCK_SG_INFOS)


def test_a_same_participant_pair_that_measures_unrelated_stays_a_conflict():
    # One participant, two SGs, but the genotypes disagree: a real sample mix-up, not the
    # two-rows-for-one-person artefact. Reclassifying it would hide the finding.
    flag = SomalierRelatednessFlag(
        category='relatedness_mismatch',
        sg_id_1='CPG014',
        sg_id_2='CPG015',
        family_external_id='FAM09',
        expected_relationship='full siblings',
        inferred_relationship='unrelated',
        relatedness=0.01,
        ibs0=9000,
        ibs2=1200,
    )

    assert not _is_same_individual(flag, MOCK_SG_INFOS)


def test_a_self_relatedness_flag_is_never_the_same_individual():
    # Both members are the same participant by definition, so if this class ever became a subclass
    # of SomalierRelatednessFlag every self-relatedness failure would be reclassified and vanish
    # from the conflicts section. The isinstance guard is what prevents that; pin it.
    flag = SomalierSelfRelatednessFlag(
        category='self_relatedness_mismatch',
        sg_id_1='CPG002',
        sg_id_2='CPG003',
        participant_external_id='PID_B',
        threshold=0.9,
        relatedness=0.61,
        ibs0=1204,
        ibs2=18337,
    )

    assert not _is_same_individual(flag, MOCK_SG_INFOS)


def test_two_sgs_with_no_resolvable_participant_are_not_the_same_individual():
    # Both participants resolve to '', which must not match each other.
    flag = SomalierRelatednessFlag(
        category='relatedness_mismatch',
        sg_id_1='CPG900',
        sg_id_2='CPG901',
        family_external_id='FAM99',
        expected_relationship='full siblings',
        inferred_relationship='identical',
        relatedness=0.998,
        ibs0=3,
        ibs2=19871,
    )

    assert not _is_same_individual(flag, {})


def test_a_sex_flag_is_never_the_same_individual():
    # Only a pairwise pedigree flag can be one, and a sex flag has no second member.
    flag = SomalierSexInferenceFlag(
        category='sex_inference_mismatch',
        provided='female',
        inferred='male',
        mean_depth=31.4,
        x_het_ratio=0.016,
        x_depth_ratio=1.02,
        y_depth_ratio=0.98,
        x_sites=4821,
        p_middling_ab=0.012,
    )

    assert not _is_same_individual(flag, MOCK_SG_INFOS)


def test_the_same_individual_section_renders_with_its_explanation():
    html = render_fixture_html()

    assert 'Same individual' in html
    # The cross-reference is the actionable part and the one claim worth pinning: it tells an
    # analyst that a genuine mismatch between two SGs of one person surfaces above as a
    # self-relatedness flag, not in this section of non-findings.
    assert 'Self-relatedness' in blurb_text(html, 'Same individual')


def test_the_same_individual_pair_is_absent_from_the_conflicts_section():
    html = render_fixture_html()

    conflicts_section = html.split('Pedigree conflicts')[1].split('Pedigree refinements')[0]
    same_individual_section = html.split('Same individual')[1].split('Resolved &mdash; past incidents')[0]

    assert 'PID_M' not in conflicts_section
    # Not merely absent from conflicts: it has to turn up in the new section. Asserting both ties
    # the extraction to the relocation, so a bug that dropped the pair entirely would still fail.
    assert 'CPG014' in same_individual_section
    assert 'FAM09' in same_individual_section


def test_the_same_individual_section_is_absent_when_there_is_nothing_in_it():
    only_refinement = [{'id': 'CPG012', 'meta': MOCK_SEQUENCING_GROUPS[11]['meta']}]
    html = render_fixture_html(only_refinement)

    assert 'Same individual' not in html


# ---------------------------------------------------------------------------
# Shared filter and sort bar
# ---------------------------------------------------------------------------
def test_every_section_is_filterable():
    html = render_fixture_html()

    # One data-section marker per rendered section, which is what the search iterates.
    for section in ('conflicts', 'refinements', 'same-individual', 'resolved'):
        assert f'data-section="{section}"' in html


def test_only_the_conflicts_section_is_category_filtered():
    # The chips are built from the conflict groups, so letting 'Sex inference' blank the whole
    # refinements section would be surprising.
    html = render_fixture_html()

    assert html.count('data-categorised="1"') == 1


def test_the_filter_bar_sits_above_every_section():
    html = render_fixture_html()

    # While the bar lived inside the conflicts section's `{% if conflict_groups %}` block, a
    # dataset with no conflicts got no search box at all. Position proves it is out of that block.
    assert html.index('id="search-input"') < html.index('Pedigree conflicts')


def test_a_dataset_of_only_refinements_still_gets_a_filter_bar():
    # The motivating case: no conflicts, but plenty to filter. The bar used to live inside the
    # conflicts block, so this dataset got no search box at all.
    many = [
        pedigree_flag(
            f'CPG{i:03d}',
            f'CPG{i + 1:03d}',
            f'FAM_R{i}',
            expected=UNSPECIFIED_RELATED,
            inferred='siblings',
            relatedness=0.4873,
            ibs0=341,
        )
        for i in range(1, MIN_GROUPS_FOR_FILTER_BAR + 3)
    ]
    groups = [{'id': f'CPG{i:03d}', 'meta': {'somalier_flags': [flag]}} for i, flag in enumerate(many, start=1)]
    infos = {
        f'CPG{i:03d}': SGInfo(
            **{
                **vars(MOCK_SG_INFOS['CPG004']),
                'sg_id': f'CPG{i:03d}',
                'family_external_id': f'FAM_R{i}',
                'participant_external_id': f'PID_R{i}',
            }
        )
        for i in range(1, MIN_GROUPS_FOR_FILTER_BAR + 4)
    }
    _, active, resolved, summary = run_pipeline(groups, infos)
    html = render_report('mock', active, resolved, summary=summary, generated_at='2026-09-14T10:00:00+00:00')

    assert summary['active_conflicts'] == 0
    assert 'All clear' in html
    assert 'id="search-input"' in html
    assert 'data-section="refinements"' in html


def test_rows_carry_the_sort_keys():
    html = render_fixture_html()

    assert 'data-newest-sg=' in html
    assert 'data-order=' in html
    assert 'Newest samples first' in html
