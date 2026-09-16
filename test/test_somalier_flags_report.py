"""
Tests for the Somalier flags report, run entirely offline against the mock fixtures.

The interesting behaviour is all in the grouping: which family a flag lands under, what happens
when a pedigree pair straddles two families, and how a family that has both active and resolved
flags gets split across the two sections.

`test_main_*` drives the command's entry point with every Metamist call and the Slack post faked,
so the wiring between those steps is covered too; the grouping tests below call the same public
functions directly because that is where the behaviour lives.
"""

from copy import deepcopy
from pathlib import Path

import pytest
from fixtures.somalier_flags import (
    CRAM_SIZE,
    CROSS_FAMILY_SUBJECT,
    FAM_CROSS,
    FAM_MIXED,
    FAM_REFINEMENT,
    FAM_SELF,
    FASTQ_SIZE,
    FIRST_SG_R1,
    FIRST_SG_R2,
    MOCK_ALL_CLEAR_SEQUENCING_GROUPS,
    MOCK_SEQUENCING_GROUPS,
    MOCK_SG_INFO_RESPONSE,
    MOCK_SG_INFOS,
    MOCK_TOTAL_SGS,
    PARTICIPANT_NO_FAMILY,
    PARTICIPANT_SELF,
    READ_DATE,
    SAME_FAMILY_SUBJECT,
    pedigree_flag,
    sex_flag,
)

from rd_qc.scripts import somalier_flags_report as report
from rd_qc.scripts.somalier_flags_report import (
    INLINE_FLAG_LIMIT,
    SgFlags,
    SGInfo,
    _extract_reads,
    _fmt_num,
    category_chips,
    collect_somalier_flags,
    flag_sg_key,
    get_previous_analysis,
    get_sg_infos,
    group_by_family,
    main,
    referenced_sg_ids,
    render_report,
    split_active_resolved,
    split_by_impact,
    summarise_flags,
    summary_message_text,
)
from rd_qc.utils import (
    UNSPECIFIED_RELATED,
    VERDICT_CONFLICT,
    VERDICT_REFINEMENT,
    SomalierRelatednessFlag,
    SomalierSelfRelatednessFlag,
    SomalierSexInferenceFlag,
)

GENERATED_AT = '2026-09-14T10:00:00+00:00'


def run_pipeline(sequencing_groups=MOCK_SEQUENCING_GROUPS, infos=MOCK_SG_INFOS):
    """The grouping pipeline, as main() runs it but without the two Metamist queries."""
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

    assert isinstance(flagged['TST001'].flags[0], SomalierSexInferenceFlag)
    assert isinstance(flagged['TST002'].flags[0], SomalierSelfRelatednessFlag)
    assert isinstance(flagged['TST004'].flags[0], SomalierRelatednessFlag)


def test_collect_skips_unrecognised_category_without_raising():
    groups = [{'id': 'TST001', 'meta': {'somalier_flags': [{'category': 'something_new', 'resolved': False}]}}]

    collected = collect_somalier_flags(groups)

    assert collected[0].flags == ()


def test_collect_skips_malformed_flag_without_raising():
    # A sex flag missing every one of its required measured fields.
    groups = [{'id': 'TST001', 'meta': {'somalier_flags': [{'category': 'sex_inference_mismatch'}]}}]

    collected = collect_somalier_flags(groups)

    assert collected[0].flags == ()


def test_collect_handles_sg_with_no_meta_at_all():
    collected = collect_somalier_flags([{'id': 'TST001', 'meta': None}])

    assert collected == [SgFlags(sg_id='TST001', flags=())]


# ---------------------------------------------------------------------------
# The SG key, and the fallback for flags recorded before it existed
# ---------------------------------------------------------------------------
def test_flag_sg_key_uses_the_recorded_key():
    # The recorded key is deliberately not what the pair IDs would derive, so this cannot pass by
    # falling through to the derivation below.
    recorded_key = 'TST900_TST901'
    raw = pedigree_flag('TST004', 'TST005', FAM_MIXED, 'parent-child', 'unrelated', verdict=VERDICT_CONFLICT)
    flag = SomalierRelatednessFlag(**{**raw, 'sequencing_group_key': recorded_key})

    assert flag_sg_key(flag, owning_sg_id='TST004') == recorded_key


def test_flag_sg_key_falls_back_to_the_pair_ids_for_legacy_flags():
    raw = pedigree_flag('TST004', 'TST005', FAM_MIXED, 'parent-child', 'unrelated', verdict=VERDICT_CONFLICT)
    flag = SomalierRelatednessFlag(**{**raw, 'sequencing_group_key': ''})

    assert flag_sg_key(flag, owning_sg_id='TST004') == 'TST004_TST005'


def test_flag_sg_key_falls_back_to_the_owning_sg_for_legacy_sex_flags():
    flag = SomalierSexInferenceFlag(**{**sex_flag('TST001', 'female', 'male'), 'sequencing_group_key': ''})

    assert flag_sg_key(flag, owning_sg_id='TST001') == 'TST001'


def test_referenced_sg_ids_includes_the_far_member_of_a_pair():
    # TST010 and TST011 only ever appear as the second member of a pair, never as an owning SG.
    flagged, _, _, _ = run_pipeline()

    referenced = referenced_sg_ids(flagged)

    assert 'TST010' in referenced
    assert 'TST011' in referenced


# ---------------------------------------------------------------------------
# Grouping
# ---------------------------------------------------------------------------
def test_cross_family_pedigree_flag_lands_in_both_families():
    _, active, _, _ = run_pipeline()
    by_label = group_by_label(active)

    fam_mixed_subjects = [f.subject for f in by_label[FAM_MIXED].flags]
    fam_cross_subjects = [f.subject for f in by_label[FAM_CROSS].flags]

    assert CROSS_FAMILY_SUBJECT in fam_mixed_subjects
    assert CROSS_FAMILY_SUBJECT in fam_cross_subjects


def test_cross_family_row_names_the_other_family_in_each_copy():
    _, active, _, _ = run_pipeline()
    by_label = group_by_label(active)

    mixed_row = next(f for f in by_label[FAM_MIXED].flags if f.subject == CROSS_FAMILY_SUBJECT)
    cross_row = next(f for f in by_label[FAM_CROSS].flags if f.subject == CROSS_FAMILY_SUBJECT)

    assert mixed_row.cross_family == FAM_CROSS
    assert cross_row.cross_family == FAM_MIXED


def test_the_cross_family_marker_reaches_the_page():
    assert 'cross-family' in render_fixture_html()


def test_a_cross_family_flag_is_rendered_twice_but_counted_once():
    _, active, _, summary = run_pipeline()

    rendered_rows = sum(len(group.flags) for group in active)
    distinct_flags = len({f.identity for g in active for f in g.flags})
    duplicated = {f.identity for g in active for f in g.flags if f.cross_family}

    # The page shows one extra row per cross-family pair, and the summary counts none of them
    # twice. Stated as a relationship so that adding any flag to the fixture cannot break it.
    assert rendered_rows == distinct_flags + len(duplicated)
    assert summary['active_flags'] == distinct_flags
    assert duplicated, 'the fixture must contain a cross-family pair for this to mean anything'


def test_same_family_pedigree_flag_lands_in_one_group_only():
    _, active, _, _ = run_pipeline()

    holders = [group.label for group in active if any(f.subject == SAME_FAMILY_SUBJECT for f in group.flags)]

    assert holders == [FAM_MIXED]


def test_a_pedigree_flag_falls_back_to_the_family_it_recorded_itself():
    # Metamist is preferred, but an SG whose family was not registered when the flag was written
    # would otherwise drop into a participant bucket and lose its family grouping entirely.
    unfamilied = SGInfo(
        sg_id='TST501',
        sg_type='genome',
        sg_technology='short-read',
        sg_platform='illumina',
        crams=[],
        fastq_pairs=[],
        other_reads=[],
        sample_external_id='EXT_501',
        sample_type='blood',
        participant_external_id='PID_501',
        family_external_id='',
    )
    flag = SomalierRelatednessFlag(
        **pedigree_flag('TST501', 'TST502', 'FAM77', 'parent-child', 'unrelated', verdict=VERDICT_CONFLICT)
    )

    groups = group_by_family([SgFlags(sg_id='TST501', flags=(flag,))], {'TST501': unfamilied})

    assert [g.label for g in groups] == ['FAM77']


def test_a_pedigree_flag_recording_no_usable_family_falls_back_to_the_participant():
    # check_pedigree writes 'unknown' when neither member had a family, which must not become a
    # family group literally labelled 'unknown'.
    unfamilied = SGInfo(
        sg_id='TST501',
        sg_type='genome',
        sg_technology='short-read',
        sg_platform='illumina',
        crams=[],
        fastq_pairs=[],
        other_reads=[],
        sample_external_id='EXT_501',
        sample_type='blood',
        participant_external_id='PID_501',
        family_external_id='',
    )
    flag = SomalierRelatednessFlag(
        **pedigree_flag('TST501', 'TST502', 'unknown', 'parent-child', 'unrelated', verdict=VERDICT_CONFLICT)
    )

    groups = group_by_family([SgFlags(sg_id='TST501', flags=(flag,))], {'TST501': unfamilied})

    keys = [group.key for group in groups]
    assert 'unknown' not in keys, "'unknown' is a placeholder, not a family"
    assert 'participant:PID_501' in keys


def test_self_relatedness_groups_under_the_family_from_metamist():
    # The flag records participant PID_B but no family; FAM03 comes from the SG info lookup.
    _, active, _, _ = run_pipeline()
    by_label = group_by_label(active)

    self_rows = [f for f in by_label[FAM_SELF].flags if f.category_key == 'self']

    # Both SGs belong to PID_B, so the subject is the one participant rather than a pair.
    assert [f.subject for f in self_rows] == [PARTICIPANT_SELF]


def test_self_relatedness_falls_back_to_a_participant_group_with_no_family():
    _, active, _, _ = run_pipeline()

    participant_groups = [group for group in active if group.key == f'participant:{PARTICIPANT_NO_FAMILY}']

    assert len(participant_groups) == 1
    assert participant_groups[0].label == f'(no family) · {PARTICIPANT_NO_FAMILY}'
    assert [f.subject for f in participant_groups[0].flags] == [PARTICIPANT_NO_FAMILY]


def test_group_carries_the_sg_info_for_both_members_of_a_pair():
    _, active, _, _ = run_pipeline()
    by_label = group_by_label(active)

    sg_ids = [info.sg_id for info in by_label[FAM_MIXED].sg_infos]

    assert set(sg_ids) >= {'TST004', 'TST005', 'TST010'}


def test_groups_sort_worst_first():
    _, active, _, _ = run_pipeline()

    totals = [sum(group.counts.values()) for group in active]

    assert len(set(totals)) > 1, 'the fixture must hold groups of differing severity'
    assert totals == sorted(totals, reverse=True)


# ---------------------------------------------------------------------------
# The per-family heading, filter chips and search text
# ---------------------------------------------------------------------------
def test_a_group_heading_names_each_category_it_holds():
    # The heading is how a reviewer decides whether a family is worth expanding.
    _, active, _, _ = run_pipeline()

    assert group_by_label(active)[FAM_SELF].count_summary == '2 flags (1 sex inference, 1 self-relatedness)'


def test_a_group_heading_with_one_flag_reads_in_the_singular():
    _, active, _, _ = run_pipeline()

    assert group_by_label(active)['FAM01'].count_summary == '1 flag (1 sex inference)'


def test_the_group_heading_reaches_the_page():
    _, active, _, _ = run_pipeline()
    heading = group_by_label(active)[FAM_MIXED].count_summary

    assert f'<div class="flag-count">{heading}</div>' in render_fixture_html()


def test_the_filter_chips_count_families_not_flags():
    # FAM02 holds two pedigree flags but is one family to review, so pedigree counts three
    # families (FAM02, FAM07, FAM08) rather than four flags.
    _, active, _, _ = run_pipeline()

    chips = category_chips(active)

    assert [(chip['key'], chip['count']) for chip in chips] == [('sex', 2), ('self', 2), ('pedigree', 3)]


def test_the_filter_chips_carry_the_collaborator_facing_labels():
    _, active, _, _ = run_pipeline()

    assert [chip['label'] for chip in category_chips(active)] == [
        'Sex inference',
        'Self-relatedness',
        'Pedigree relatedness',
    ]


def test_a_category_with_no_families_gets_no_chip():
    _, active, _, _ = run_pipeline(MOCK_ALL_CLEAR_SEQUENCING_GROUPS)

    assert category_chips(active) == []


@pytest.mark.parametrize(
    'term',
    [
        'fam02',  # the family label
        'fam07',  # the far family of the cross-family pair
        'pid_c',  # a participant external ID
        'tst004',  # an SG ID
        'ext_d',  # a sample external ID
        'parent-child',  # the recorded relationship
    ],
)
def test_a_family_is_findable_by_any_identifier_it_involves(term):
    # The filter box matches against this blob, so anything a collaborator might paste has to be
    # in it — including the identifiers of the far family of a cross-family pair.
    _, active, _, _ = run_pipeline()

    assert term in group_by_label(active)[FAM_MIXED].search_blob


def test_the_search_text_reaches_the_page():
    _, active, _, _ = run_pipeline()
    blob = group_by_label(active)[FAM_MIXED].search_blob

    assert f'data-search="{blob}"' in render_fixture_html()


# ---------------------------------------------------------------------------
# Active / resolved split
# ---------------------------------------------------------------------------
def test_family_with_both_kinds_appears_in_both_sections():
    _, active, resolved, _ = run_pipeline()

    assert FAM_MIXED in group_by_label(active)
    assert FAM_MIXED in group_by_label(resolved)


def test_each_section_shows_only_its_own_flags():
    _, active, resolved, _ = run_pipeline()

    assert active and resolved, 'the fixture must supply both an active and a resolved flag'
    assert all(not f.resolved for group in active for f in group.flags)
    assert all(f.resolved for group in resolved for f in group.flags)


def test_counts_are_recomputed_per_section():
    _, active, resolved, _ = run_pipeline()

    # FAM02 has two active pedigree flags, and one resolved pedigree plus one resolved sex.
    assert group_by_label(active)[FAM_MIXED].counts == {'pedigree': 2}
    assert group_by_label(resolved)[FAM_MIXED].counts == {'sex': 1, 'pedigree': 1}


# ---------------------------------------------------------------------------
# Summary counts
# ---------------------------------------------------------------------------
def test_summary_counts():
    flagged, _, _, _ = run_pipeline()

    summary = summarise_flags(flagged, total_sgs=MOCK_TOTAL_SGS, families_affected=6)

    assert summary == {
        'total_sgs': 13,
        'active_flags': 7,
        'active_by_category': {'sex': 2, 'self': 2, 'pedigree': 3},
        # Six of the seven are real disagreements; FAM08's unspecified/siblings pair is not.
        # That one is a pedigree flag, so pedigree drops to 2 once refinements come out.
        'active_conflicts_by_category': {'sex': 2, 'self': 2, 'pedigree': 2},
        'active_conflicts': 6,
        'active_refinements': 1,
        'families_affected': 6,
        'resolved_flags': 2,
    }


def test_the_families_affected_count_is_the_one_the_caller_supplied():
    # summarise_flags counts flags, not groups, so the group count is passed in. A value that is
    # not the fixture's real group count proves it is carried rather than recomputed.
    flagged, _, _, _ = run_pipeline()

    summary = summarise_flags(flagged, total_sgs=MOCK_TOTAL_SGS, families_affected=41)

    assert summary['families_affected'] == 41


def test_a_flag_recorded_before_verdict_existed_counts_as_a_conflict():
    # Records written before the field was added carry ''. Treating those as refinements would
    # silently move every pre-existing pedigree flag into the de-emphasised section.
    legacy_flag = pedigree_flag('TST004', 'TST005', FAM_MIXED, 'parent-child', 'unrelated', verdict='')
    flagged = list(collect_somalier_flags([{'id': 'TST004', 'meta': {'somalier_flags': [legacy_flag]}}]))

    summary = summarise_flags(flagged, total_sgs=1, families_affected=1)

    assert (summary['active_conflicts'], summary['active_refinements']) == (1, 0)


def test_a_flag_recorded_before_verdict_existed_renders_as_a_conflict():
    legacy_flag = pedigree_flag('TST004', 'TST005', FAM_MIXED, 'parent-child', 'unrelated', verdict='')
    flagged = list(collect_somalier_flags([{'id': 'TST004', 'meta': {'somalier_flags': [legacy_flag]}}]))
    active, _ = split_active_resolved(group_by_family(flagged, MOCK_SG_INFOS), MOCK_SG_INFOS)

    conflicts, refinements = split_by_impact(active, MOCK_SG_INFOS)

    assert [g.label for g in conflicts] == [FAM_MIXED]
    assert refinements == []


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
        (1.5, '1.5'),
        (0.9, '0.9'),
        (2.0, '2'),
        (None, '—'),
        ('', '—'),
        ('unknown', 'unknown'),
        # A bool is an int in Python, so without the explicit guard True would format as '1'.
        (True, 'True'),
    ],
)
def test_fmt_num(value, expected):
    assert _fmt_num(value) == expected


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------
def render_fixture_html(sequencing_groups=MOCK_SEQUENCING_GROUPS) -> str:
    _, active, resolved, summary = run_pipeline(sequencing_groups)
    return render_report('mock-dataset', active, resolved, summary=summary, generated_at=GENERATED_AT)


def test_render_emits_a_row_for_every_flag_in_every_group():
    # A positive floor on the rendered page: a template that silently stopped emitting rows, or
    # rendered an empty document, fails here where an "absence of leftovers" check would not.
    _, active, resolved, _ = run_pipeline()
    expected_rows = sum(len(g.flags) for g in active) + sum(len(g.flags) for g in resolved)

    html = render_fixture_html()

    assert html.count('<tr') >= expected_rows
    assert '{{' not in html, 'unrendered template syntax reached the page'
    assert '{%' not in html, 'unrendered template syntax reached the page'


@pytest.mark.parametrize(
    'result_line',
    [
        'expected parent-child / measured unrelated',
        'provided unknown / inferred male',
        'relatedness 0.62 (expected ~1.0, threshold 0.9)',
    ],
)
def test_each_category_renders_its_own_result_line(result_line):
    assert result_line in render_fixture_html()


def test_the_resolved_section_appears_when_something_is_resolved():
    assert 'Resolved &mdash; past incidents' in render_fixture_html()


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


# ---------------------------------------------------------------------------
# Capping the inline flag lines
# ---------------------------------------------------------------------------
# One family, more flags than it will show inline. A real family can carry 50+ pedigree
# mismatches, which rendered inline would bury every other family on the page.
CAPPED_FAMILY = 'FAM09'
FLAG_COUNT = INLINE_FLAG_LIMIT + 3


def capped_family_inputs() -> tuple[list, dict]:
    """One over-full family: FLAG_COUNT pedigree flags between FLAG_COUNT + 1 distinct SGs."""
    flags = [
        pedigree_flag(
            f'TST{i:03d}',
            f'TST{i + 1:03d}',
            CAPPED_FAMILY,
            expected='unrelated',
            inferred='full siblings',
            verdict=VERDICT_CONFLICT,
        )
        for i in range(1, FLAG_COUNT + 1)
    ]
    infos = {
        f'TST{i:03d}': SGInfo(
            sg_id=f'TST{i:03d}',
            sg_type='genome',
            sg_technology='short-read',
            sg_platform='illumina',
            crams=[],
            fastq_pairs=[],
            other_reads=[],
            sample_external_id=f'EXT_{i:03d}',
            sample_type='blood',
            # Distinct per SG, so a row that names the wrong pair is visible.
            participant_external_id=f'PID_{i:03d}',
            family_external_id=CAPPED_FAMILY,
        )
        for i in range(1, FLAG_COUNT + 2)
    }
    return flags, infos


def capped_family_render() -> tuple[list, str]:
    flags, infos = capped_family_inputs()
    flagged = [sf for sf in collect_somalier_flags([{'id': 'TST001', 'meta': {'somalier_flags': flags}}]) if sf.flags]
    active, resolved = split_active_resolved(group_by_family(flagged, infos), infos)
    summary = summarise_flags(flagged, total_sgs=1, families_affected=len(active))
    return active, render_report('mock', active, resolved, summary=summary, generated_at=GENERATED_AT)


def test_a_groups_total_counts_every_flag_not_just_the_inline_ones():
    active, _ = capped_family_render()

    assert active[0].total == FLAG_COUNT


def test_the_inline_list_is_capped_with_a_more_link():
    _, html = capped_family_render()

    assert f'+{FLAG_COUNT - INLINE_FLAG_LIMIT} more &mdash; click to expand' in html


def test_every_capped_flag_still_reaches_the_detail_table():
    # The inline list stops at the cap, but the expandable detail table below it carries a row per
    # flag, naming both members. Without that, the capped flags would be invisible on the page.
    _, html = capped_family_render()

    involved_sgs = [f'TST{i:03d}' for i in range(1, FLAG_COUNT + 2)]

    assert [sg for sg in involved_sgs if f'<div class="member-sgid">{sg}</div>' not in html] == []


# ---------------------------------------------------------------------------
# Conflicts and refinements
# ---------------------------------------------------------------------------
def test_split_by_impact_separates_conflicts_from_refinements():
    _, active, _, _ = run_pipeline()

    conflicts, refinements = split_by_impact(active, {info.sg_id: info for g in active for info in g.sg_infos})

    assert conflicts and refinements, 'the fixture must supply both a conflict and a refinement'
    # FAM08's only flag records no path between the pair, so it is solely a refinement.
    assert FAM_REFINEMENT not in group_by_label(conflicts)
    assert FAM_REFINEMENT in group_by_label(refinements)
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
    only_refinement = [
        {
            'id': 'TST012',
            'meta': {
                'somalier_flags': [
                    pedigree_flag(
                        'TST012',
                        'TST013',
                        FAM_REFINEMENT,
                        expected=UNSPECIFIED_RELATED,
                        inferred='siblings',
                        verdict=VERDICT_REFINEMENT,
                    ),
                ],
            },
        },
    ]

    html = render_fixture_html(only_refinement)

    assert 'All clear' in html
    assert 'Pedigree refinements' in html


# ---------------------------------------------------------------------------
# Slack summary message
# ---------------------------------------------------------------------------
FLAGS_URL = 'https://main-web.populationgenomics.org.au/mock/flags.html'
SOMALIER_URL = 'https://main-web.populationgenomics.org.au/mock/somalier.html'
PREVIOUS_REPORT_DATE = '2026-08-20'


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


def previous_report(summary: dict, completed: str = f'{PREVIOUS_REPORT_DATE}T03:00:00+00:00') -> dict:
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
    assert 'active conflict' not in text


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
    # Spelled out rather than derived from `summary`, so both sides of every delta are in view.
    previous = previous_report(
        {
            'total_sgs': 11,
            'active_conflicts': 4,
            'active_refinements': 0,
            'families_affected': 5,
            'resolved_flags': 0,
        }
    )

    text = message_for(summary, previous)

    assert f'Changes since the last report on {PREVIOUS_REPORT_DATE}' in text
    assert '+2 new sequencing groups' in text
    assert '+1 additional family flagged' in text
    assert '+2 new conflicts' in text
    assert '+1 new pedigree refinement' in text


def test_message_reports_no_change_when_the_summary_is_identical():
    _, _, _, summary = run_pipeline()

    text = message_for(summary, previous_report(dict(summary)))

    assert f'No change since the last report on {PREVIOUS_REPORT_DATE}' in text
    assert 'Changes since' not in text


def test_message_reports_flags_that_have_been_fixed():
    _, _, _, summary = run_pipeline()
    previous = previous_report({**summary, 'active_conflicts': 9, 'resolved_flags': 0})

    text = message_for(summary, previous)

    assert '3 fewer conflicts' in text
    assert '2 more flags resolved' in text


def test_message_survives_a_previous_summary_from_before_the_conflict_split():
    _, _, _, summary = run_pipeline()
    # Reports registered before conflicts/refinements existed only carry these keys. The missing
    # ones have to read as zero rather than raising, which makes this the change branch.
    previous = previous_report({'total_sgs': 13, 'active_flags': 7, 'resolved_flags': 2})

    text = message_for(summary, previous)

    assert f'Changes since the last report on {PREVIOUS_REPORT_DATE}' in text
    assert '+6 new conflicts' in text
    assert '+1 new pedigree refinement' in text


# ---------------------------------------------------------------------------
# Wording of an unspecified expectation
# ---------------------------------------------------------------------------
def refinement_row():
    """FAM08's only flag: no recorded path between the pair, measured as siblings."""
    _, active, _, _ = run_pipeline()
    flags = group_by_label(active)[FAM_REFINEMENT].flags
    assert len(flags) == 1, 'FAM08 must hold exactly one flag for this row to be unambiguous'
    return flags[0]


def test_an_unspecified_expectation_reads_as_no_relationship_provided():
    # 'expected related at unknown level / measured siblings' is accurate but unclear: the point is
    # that the pedigree has no data here, not that it asserted something vague.
    assert refinement_row().result == 'No relationship provided / measured siblings'


def test_the_detail_table_uses_the_same_wording():
    assert refinement_row().expected == 'No relationship provided'


def test_a_stated_expectation_is_left_alone():
    _, active, _, _ = run_pipeline()
    row = group_by_label(active)[FAM_MIXED].flags[0]

    assert row.result == 'expected parent-child / measured unrelated'
    assert row.expected == 'parent-child'


def test_the_raw_relationship_is_still_searchable():
    # A collaborator pasting the peddy string, or a saved filter, must still find the row. The
    # stored value is peddy's vocabulary and part of the flag's reconciliation identity, so only
    # the display may be relabelled.
    assert UNSPECIFIED_RELATED in refinement_row().search_blob


def test_the_refinements_blurb_describes_the_case_that_actually_occurs():
    html = render_fixture_html()
    blurb = html.split('Pedigree refinements', 1)[1].split('</section>', 1)[0]

    # 'siblings' -> 'full siblings' is satisfied by EXPECTED_DEGREES now, so it never reaches the
    # refinements section and must not be described as the common case.
    assert 'whenever a parent is missing' in blurb
    assert 'full siblings' not in blurb


# ---------------------------------------------------------------------------
# Glance line leads with the identifiers collaborators use
# ---------------------------------------------------------------------------
def test_a_pedigree_pair_leads_with_participants_and_demotes_the_sg_ids():
    _, active, _, _ = run_pipeline()
    row = next(f for f in group_by_label(active)[FAM_MIXED].flags if f.category_key == 'pedigree')

    assert row.subject == SAME_FAMILY_SUBJECT
    assert row.subject_detail == 'TST004 ↔ TST005'


def test_a_self_relatedness_pair_leads_with_the_one_participant():
    _, active, _, _ = run_pipeline()
    row = next(f for f in group_by_label(active)[FAM_SELF].flags if f.category_key == 'self')

    assert row.subject == PARTICIPANT_SELF
    assert row.subject_detail == 'TST002 ↔ TST003'


def test_a_pair_with_no_known_participants_falls_back_to_the_sg_ids():
    # Nothing to promote, so the SG ids stay in the lead rather than leaving the line blank.
    raw = pedigree_flag('TST404', 'TST405', 'FAM99', 'parent-child', 'unrelated', verdict=VERDICT_CONFLICT)
    flag = SomalierRelatednessFlag(**raw)
    rows = group_by_family([SgFlags(sg_id='TST404', flags=(flag,))], {})

    assert rows[0].flags[0].subject == 'TST404 ↔ TST405'
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
    assert f'{FIRST_SG_R1}</span>' in html
    assert '&nbsp;/&nbsp;' not in html


def test_the_rendered_read_lines_show_size_and_date_when_known():
    html = render_fixture_html()

    # Unlabelled: a GiB figure and a date read as themselves.
    assert f'<span class="read-meta">{FASTQ_SIZE}</span>' in html
    assert f'<span class="read-meta">{READ_DATE}</span>' in html


def test_a_read_with_no_size_or_date_renders_the_name_alone():
    # The fixture's R2 has neither, so it must not emit an empty meta span.
    html = render_fixture_html()
    r2_line = next(line for line in html.splitlines() if FIRST_SG_R2 in line)

    assert 'read-meta' not in r2_line


# ---------------------------------------------------------------------------
# Parsing Metamist's SG info response
# ---------------------------------------------------------------------------
def test_get_sg_infos_builds_the_same_infos_the_fixture_mirrors(monkeypatch):
    """
    The gate on MOCK_SG_INFOS: it claims to be what `get_sg_infos` returns, so prove it.

    Without this, every render test could pass against a hand-built mirror while the real parser
    mis-read Metamist's externalIds-keyed-by-'' quirk, its per-assay read grouping or its byte
    sizes, and nothing would fail.
    """
    monkeypatch.setattr(report, 'query', lambda _q, **_kwargs: MOCK_SG_INFO_RESPONSE)

    infos = get_sg_infos(list(MOCK_SG_INFOS))

    assert infos == MOCK_SG_INFOS


def test_get_sg_infos_reads_the_primary_external_ids_and_family(monkeypatch):
    # Named separately from the whole-dict comparison above so a failure says which field moved.
    monkeypatch.setattr(report, 'query', lambda _q, **_kwargs: MOCK_SG_INFO_RESPONSE)

    info = get_sg_infos(['TST001'])['TST001']

    assert (info.sample_external_id, info.participant_external_id, info.family_external_id) == (
        'EXT_A',
        'PID_A',
        'FAM01',
    )


def test_get_sg_infos_formats_the_recorded_byte_size(monkeypatch):
    monkeypatch.setattr(report, 'query', lambda _q, **_kwargs: MOCK_SG_INFO_RESPONSE)

    info = get_sg_infos(['TST001'])['TST001']

    assert (info.crams[0].size, info.fastq_pairs[0][0].size) == (CRAM_SIZE, FASTQ_SIZE)


def test_get_sg_infos_leaves_an_sg_with_no_family_unfamilied(monkeypatch):
    # Metamist returns families: [], which must read as '' rather than raising on families[0].
    monkeypatch.setattr(report, 'query', lambda _q, **_kwargs: MOCK_SG_INFO_RESPONSE)

    assert get_sg_infos(['TST009'])['TST009'].family_external_id == ''


def test_get_sg_infos_skips_the_query_entirely_for_an_empty_sg_list(monkeypatch):
    def fail_on_contact(_q, **_kwargs: object) -> dict:
        raise AssertionError('must not query Metamist for an empty SG list')

    monkeypatch.setattr(report, 'query', fail_on_contact)

    assert get_sg_infos([]) == {}


# ---------------------------------------------------------------------------
# Finding the previous report to compare against
# ---------------------------------------------------------------------------
def analyses_response(analyses: list[dict]) -> dict:
    return {'project': {'analyses': analyses}}


def analysis(analysis_id: int, completed: str, summary: dict | None) -> dict:
    return {
        'id': analysis_id,
        'outputs': f'gs://bucket/report-{analysis_id}.html',
        'timestampCompleted': completed,
        'meta': {'summary': summary} if summary is not None else {},
    }


def test_a_single_existing_analysis_is_not_a_previous_report(monkeypatch):
    # The only analysis on record is the one this run just registered.
    monkeypatch.setattr(
        report,
        'query',
        lambda _q, **_kwargs: analyses_response([analysis(1, '2026-09-14T10:00:00+00:00', {'active_flags': 7})]),
    )

    assert get_previous_analysis('mock-dataset', {'stage': 'GenerateSomalierFlagsReport'}) is None


def test_no_existing_analyses_means_no_previous_report(monkeypatch):
    monkeypatch.setattr(report, 'query', lambda _q, **_kwargs: analyses_response([]))

    assert get_previous_analysis('mock-dataset', {'stage': 'GenerateSomalierFlagsReport'}) is None


def test_the_second_most_recent_analysis_is_the_previous_report(monkeypatch):
    # Deliberately out of order in the response, so the sort is what picks the right one.
    monkeypatch.setattr(
        report,
        'query',
        lambda _q, **_kwargs: analyses_response(
            [
                analysis(2, '2026-08-20T03:00:00+00:00', {'active_flags': 5}),
                analysis(3, '2026-09-14T10:00:00+00:00', {'active_flags': 7}),
                analysis(1, '2026-07-01T03:00:00+00:00', {'active_flags': 3}),
            ]
        ),
    )

    previous = get_previous_analysis('mock-dataset', {'stage': 'GenerateSomalierFlagsReport'})

    assert previous is not None, 'two prior analyses exist, so one of them is the previous report'
    assert previous['id'] == 2


def test_a_previous_analysis_without_a_summary_is_ignored(monkeypatch):
    # Nothing to diff against, so the message must omit the change section rather than compare
    # against zeroes and claim everything is new.
    monkeypatch.setattr(
        report,
        'query',
        lambda _q, **_kwargs: analyses_response(
            [
                analysis(2, '2026-08-20T03:00:00+00:00', summary=None),
                analysis(3, '2026-09-14T10:00:00+00:00', {'active_flags': 7}),
            ]
        ),
    )

    assert get_previous_analysis('mock-dataset', {'stage': 'GenerateSomalierFlagsReport'}) is None


# ---------------------------------------------------------------------------
# main(), with Metamist and Slack faked
# ---------------------------------------------------------------------------
DATASET = 'mock-dataset'
RESOLVED_DATASET = 'mock-dataset-test'
SEQ_TYPE = 'genome'
SEQ_TECH = 'short-read'


@pytest.fixture
def report_run(monkeypatch, tmp_path):
    """
    Everything `main` reaches outside the process, faked, plus a record of what it did.

    The two query responses are the recorded-shape fixtures, so `main` runs the real
    `collect_somalier_flags`/`get_sg_infos`/grouping/render path end to end.
    """
    run: dict = {'queries': [], 'registered': None, 'slack': None}

    run['config'] = {
        ('workflow', 'sequencing_type'): SEQ_TYPE,
        ('workflow', 'sequencing_technology'): SEQ_TECH,
        ('somalier_flags_report', 'send_to_slack'): True,
    }

    def fake_config_retrieve(key, **_kwargs: object) -> object:
        """
        Resolves the keys main() reads, and raises on anything else rather than guessing.

        Production defaults are ignored deliberately: every key main() reads is spelled out above,
        so a new one shows up as a failure here rather than silently taking its default.
        """
        if tuple(key) not in run['config']:
            raise AssertionError(f'unexpected config key {key!r}')
        return run['config'][tuple(key)]

    def fake_query(query_obj, variables=None) -> dict:
        run['queries'].append((query_obj, variables))
        if query_obj is report.DATASET_SGS_QUERY:
            return {'project': {'sequencingGroups': MOCK_SEQUENCING_GROUPS}}
        if query_obj is report.SGS_INFO_QUERY:
            return MOCK_SG_INFO_RESPONSE
        if query_obj is report.EXISTING_ANALYSES_QUERY:
            return analyses_response([])
        raise AssertionError('main() issued an unrecognised query')

    monkeypatch.setattr(report, 'query', fake_query)
    monkeypatch.setattr(report, 'config_retrieve', fake_config_retrieve)
    monkeypatch.setattr(report, 'dataset_for_access_level', lambda _name: RESOLVED_DATASET)

    def fake_create_new(**kwargs: object) -> None:
        # The real client serialises the payload during the call. main() then pops 'summary' from
        # that same dict, so capturing by reference here would not see what was actually sent.
        run['registered'] = deepcopy(kwargs)

    monkeypatch.setattr(report, 'create_new', fake_create_new)
    monkeypatch.setattr(report, 'send_message', lambda text: run.__setitem__('slack', text))

    run['output_html'] = str(tmp_path / 'report-2026-09-14.html')
    run['base_output_html'] = str(tmp_path / 'report.html')
    return run


def run_main(report_run) -> dict:
    main(
        dataset=DATASET,
        output_html=report_run['output_html'],
        base_output_html=report_run['base_output_html'],
        flags_html_url=FLAGS_URL,
        somalier_html_url=SOMALIER_URL,
    )
    return report_run


def test_main_writes_the_report_to_both_the_timestamped_and_fixed_paths(report_run, tmp_path):
    run = run_main(report_run)

    timestamped = (tmp_path / 'report-2026-09-14.html').read_text()
    fixed = (tmp_path / 'report.html').read_text()

    assert 'Pedigree conflicts' in fixed, 'the fixed path must hold a rendered report'
    assert timestamped == fixed
    assert run['registered'] is not None


def test_main_renders_the_flags_it_read_from_metamist(report_run):
    run_main(report_run)

    html = Path(report_run['output_html']).read_text()

    # The cross-family pair and the refinement are the two grouping branches main() has to wire
    # the SG-info query into; both need participant IDs that only that query supplies.
    assert CROSS_FAMILY_SUBJECT in html
    assert FAM_REFINEMENT in html


def test_main_registers_the_analysis_against_every_dataset_sg(report_run):
    run = run_main(report_run)

    registered = run['registered']
    assert registered['project'] == RESOLVED_DATASET
    assert registered['output'] == report_run['output_html']
    assert registered['analysis_type'] == 'web'
    assert set(registered['sgs']) == {sg['id'] for sg in MOCK_SEQUENCING_GROUPS}


def test_main_records_the_summary_in_the_analysis_meta(report_run):
    # The next run diffs against this, so the counts have to be the ones that were rendered.
    run = run_main(report_run)

    summary = run['registered']['meta']['summary']

    assert summary['total_sgs'] == MOCK_TOTAL_SGS
    assert summary['active_conflicts'] == 6
    assert summary['active_refinements'] == 1
    assert summary['families_affected'] == 6


def test_main_looks_for_the_previous_report_without_the_summary_in_the_filter(report_run):
    # The meta filter has to match a previous analysis, whose summary differs by construction.
    # Leaving 'summary' in the filter would mean no previous report is ever found.
    run = run_main(report_run)

    meta_filters = [v['metaFilter'] for q, v in run['queries'] if q is report.EXISTING_ANALYSES_QUERY]

    assert len(meta_filters) == 1, 'main() must look for a previous report exactly once'
    assert 'summary' not in meta_filters[0]
    assert meta_filters[0]['stage'] == report.STAGE_NAME
    assert meta_filters[0]['dataset'] == RESOLVED_DATASET


def test_main_queries_metamist_with_the_access_level_resolved_dataset(report_run):
    run = run_main(report_run)

    dataset_query = next(v for q, v in run['queries'] if q is report.DATASET_SGS_QUERY)

    assert dataset_query['dataset'] == RESOLVED_DATASET
    assert (dataset_query['seqType'], dataset_query['seqTech']) == (SEQ_TYPE, SEQ_TECH)


def test_main_posts_the_summary_to_slack(report_run):
    run = run_main(report_run)

    assert run['slack'] is not None, 'main() must post the summary when send_to_slack is on'
    # The resolved name, so the message says which project was actually read.
    assert run['slack'].startswith(f'*[{RESOLVED_DATASET}]* <{FLAGS_URL}|')
    assert '*6 active conflicts*' in run['slack']


def test_main_does_not_post_to_slack_when_the_config_turns_it_off(report_run):
    report_run['config'][('somalier_flags_report', 'send_to_slack')] = False

    run = run_main(report_run)

    assert run['slack'] is None
    assert run['registered'] is not None, 'the analysis is still registered with Slack off'
