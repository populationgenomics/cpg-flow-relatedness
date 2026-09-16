"""
Tests for the per-participant self-relatedness check.

Two SGs from the same person must measure as effectively the same genome. When they do not, one of
them is mislabelled, so this check is the one that catches a sample swap within a participant.

The three functions covered here are pure: a somalier pairs.tsv in, and the Slack alert plus the
flags to record out. `run` itself only adds Metamist registration and file copying on top.
"""

import pytest
from fixtures.somalier_tsv import pair_row, write_pairs

from rd_qc.scripts.check_self_relatedness import (
    build_alert_and_flags,
    find_existing_analysis,
    read_low_relatedness_pairs,
)

# Two SGs of one participant are expected at ~1.0, so anything under this is a mismatch.
THRESHOLD = 0.9
DATASET = 'test-dataset'
PARTICIPANT = 'PID_X'
# Distinct from the external ID, since older analyses key on one and newer ones on the other.
PARTICIPANT_ID = 4242
HTML_URL = 'https://main-web.populationgenomics.org.au/test/somalier.html'


# ---------------------------------------------------------------------------
# Reading the pairs file
# ---------------------------------------------------------------------------
def test_only_pairs_below_the_threshold_are_returned(tmp_path):
    pairs_path = write_pairs(
        tmp_path,
        [
            pair_row('TST001', 'TST002', relatedness=0.42),
            pair_row('TST001', 'TST003', relatedness=0.98),
        ],
    )

    low = read_low_relatedness_pairs(pairs_path, THRESHOLD)

    assert [(p['sample_a'], p['sample_b']) for p in low] == [('TST001', 'TST002')]


def test_a_pair_exactly_at_the_threshold_is_not_low(tmp_path):
    # The comparison is strict, so a pair sitting on the threshold passes.
    pairs_path = write_pairs(tmp_path, [pair_row('TST001', 'TST002', relatedness=THRESHOLD)])

    assert read_low_relatedness_pairs(pairs_path, THRESHOLD) == []


def test_the_measured_values_come_through_with_the_pair(tmp_path):
    pairs_path = write_pairs(tmp_path, [pair_row('TST001', 'TST002', relatedness=0.42, ibs0=1204)])

    low = read_low_relatedness_pairs(pairs_path, THRESHOLD)

    assert low[0]['relatedness'] == 0.42
    assert int(low[0]['ibs0']) == 1204


def test_a_pairs_file_with_no_rows_reads_as_empty(tmp_path):
    pairs_path = write_pairs(tmp_path, [])

    assert read_low_relatedness_pairs(pairs_path, THRESHOLD) == []


def test_a_missing_pairs_file_reads_as_none_rather_than_empty(tmp_path):
    # `run` distinguishes the two: None means somalier produced nothing and the check exits early,
    # while [] means it ran and found no mismatch, which registers a passing analysis.
    missing = str(tmp_path / 'never-written.pairs.tsv')

    assert read_low_relatedness_pairs(missing, THRESHOLD) is None


# ---------------------------------------------------------------------------
# Finding a participant's existing analysis
# ---------------------------------------------------------------------------
def analysis(meta: dict) -> dict:
    return {'id': 501, 'meta': meta}


@pytest.mark.parametrize(
    'meta',
    [
        # Oldest records key on the internal integer ID.
        {'participant_id': PARTICIPANT_ID},
        # Then the external ID moved under the same key.
        {'participant_id': PARTICIPANT},
        # Newest records use a key of their own.
        {'participant_external_id': PARTICIPANT},
    ],
)
def test_an_existing_analysis_is_found_under_any_recorded_participant_key(meta):
    analyses = [analysis({'participant_id': 9999}), analysis(meta)]

    found = find_existing_analysis(analyses, PARTICIPANT_ID, PARTICIPANT)

    assert found is not None, 'the participant is recorded, so its analysis must be found'
    assert found['meta'] == meta


def test_a_participant_with_no_recorded_analysis_is_not_found():
    # A different participant's records must not be claimed as this one's, or a passing check
    # would silently skip registering its own analysis.
    analyses = [analysis({'participant_id': 9999}), analysis({'participant_external_id': 'PID_OTHER'})]

    assert find_existing_analysis(analyses, PARTICIPANT_ID, PARTICIPANT) is None


def test_no_analyses_at_all_is_not_found():
    assert find_existing_analysis([], PARTICIPANT_ID, PARTICIPANT) is None


# ---------------------------------------------------------------------------
# Building the alert and the flags
# ---------------------------------------------------------------------------
def low_pair(sample_a: str, sample_b: str, relatedness: float = 0.42, ibs0: int = 1204, ibs2: int = 18337) -> dict:
    """One entry in the shape `read_low_relatedness_pairs` yields, with ibs0/ibs2 as read strings."""
    return {
        'sample_a': sample_a,
        'sample_b': sample_b,
        'relatedness': relatedness,
        'ibs0': str(ibs0),
        'ibs2': str(ibs2),
    }


def test_a_flag_is_recorded_against_the_first_sg_of_the_pair_only():
    # The flags report dedupes pairwise flags on a sorted key, so recording against both members
    # would render the same mismatch twice.
    _, flags = build_alert_and_flags(DATASET, PARTICIPANT, THRESHOLD, HTML_URL, [low_pair('TST001', 'TST002')])

    assert list(flags) == ['TST001']


def test_a_pair_reported_backwards_is_sorted_before_it_is_recorded():
    _, flags = build_alert_and_flags(DATASET, PARTICIPANT, THRESHOLD, HTML_URL, [low_pair('TST002', 'TST001')])

    flag = flags['TST001'][0]
    assert (flag.sg_id_1, flag.sg_id_2) == ('TST001', 'TST002')


def test_the_flag_carries_the_measured_values_as_numbers():
    # ibs0/ibs2 arrive as strings from the TSV reader, and Metamist meta has to hold numbers.
    _, flags = build_alert_and_flags(
        DATASET,
        PARTICIPANT,
        THRESHOLD,
        HTML_URL,
        [low_pair('TST001', 'TST002', relatedness=0.42, ibs0=1204, ibs2=18337)],
    )

    flag = flags['TST001'][0]
    assert (flag.relatedness, flag.ibs0, flag.ibs2) == (0.42, 1204, 18337)
    assert (flag.threshold, flag.participant_external_id) == (THRESHOLD, PARTICIPANT)


def test_the_alert_names_the_participant_and_the_threshold():
    text, _ = build_alert_and_flags(DATASET, PARTICIPANT, THRESHOLD, HTML_URL, [low_pair('TST001', 'TST002')])

    assert f'*[{DATASET}]*' in text
    assert f'Self-relatedness check failed for participant {PARTICIPANT}' in text
    assert f'threshold: {THRESHOLD}' in text


def test_the_alert_lists_every_failing_pair_with_its_measurements():
    text, _ = build_alert_and_flags(
        DATASET,
        PARTICIPANT,
        THRESHOLD,
        HTML_URL,
        [low_pair('TST001', 'TST002', relatedness=0.42), low_pair('TST001', 'TST003', relatedness=0.11)],
    )

    assert 'TST001 - TST002: relatedness=0.42' in text
    assert 'TST001 - TST003: relatedness=0.11' in text


def test_the_alert_links_the_report_when_there_is_a_url():
    text, _ = build_alert_and_flags(DATASET, PARTICIPANT, THRESHOLD, HTML_URL, [low_pair('TST001', 'TST002')])

    assert f'<{HTML_URL}|Self-relatedness check failed for participant {PARTICIPANT}>' in text


def test_the_alert_reads_as_plain_text_when_there_is_no_url():
    # The check still runs when no web report was produced, and the alert must not emit an
    # empty Slack link in that case.
    text, _ = build_alert_and_flags(DATASET, PARTICIPANT, THRESHOLD, '', [low_pair('TST001', 'TST002')])

    assert '<|' not in text
    assert f'*[{DATASET}]* Self-relatedness check failed for participant {PARTICIPANT}' in text


def test_a_pairs_file_feeds_straight_into_the_alert(tmp_path):
    # The two halves share the pair dict's shape, including ibs0/ibs2 arriving as strings, so
    # they are exercised together rather than only against hand-built dicts.
    pairs_path = write_pairs(tmp_path, [pair_row('TST002', 'TST001', relatedness=0.42, ibs0=1204)])

    low = read_low_relatedness_pairs(pairs_path, THRESHOLD)
    text, flags = build_alert_and_flags(DATASET, PARTICIPANT, THRESHOLD, HTML_URL, low)

    assert 'TST001 - TST002: relatedness=0.42, ibs0=1204' in text
    assert flags['TST001'][0].ibs0 == 1204
