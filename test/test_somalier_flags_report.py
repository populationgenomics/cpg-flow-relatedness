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
    _fmt_num,
    collect_somalier_flags,
    flag_sg_key,
    group_by_family,
    referenced_sg_ids,
    render_report,
    split_active_resolved,
    summarise_flags,
)
from rd_qc.utils import SomalierRelatednessFlag, SomalierSelfRelatednessFlag, SomalierSexInferenceFlag

CROSS_FAMILY_SUBJECT = 'CPG004 ↔ CPG010'


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

    # Seven rows on the page, but only six distinct flags: the cross-family pair is shown twice.
    assert rendered_rows == 7
    assert summary['active_flags'] == 6


def test_same_family_pedigree_flag_lands_in_one_group_only():
    _, active, _, _ = run_pipeline()
    same_family_subject = 'CPG004 ↔ CPG005'

    holders = [group.label for group in active if any(f.subject == same_family_subject for f in group.flags)]

    assert holders == ['FAM02']


def test_self_relatedness_groups_under_the_family_from_metamist():
    # The flag records participant PID_B but no family; FAM03 comes from the SG info lookup.
    _, active, _, _ = run_pipeline()
    by_label = group_by_label(active)

    self_rows = [f for f in by_label['FAM03'].flags if f.category_key == 'self']

    assert [f.subject for f in self_rows] == ['CPG002 ↔ CPG003']


def test_self_relatedness_falls_back_to_a_participant_group_with_no_family():
    _, active, _, _ = run_pipeline()

    participant_groups = [group for group in active if group.key == 'participant:PID_G']

    assert len(participant_groups) == 1
    assert participant_groups[0].label == '(no family) · PID_G'
    assert [f.subject for f in participant_groups[0].flags] == ['CPG009 ↔ CPG011']


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
        'total_sgs': 11,
        'active_flags': 6,
        'active_by_category': {'sex': 2, 'self': 2, 'pedigree': 2},
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
        'expected parent-child / inferred unrelated',
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
    assert 'Families with flags' in html


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
