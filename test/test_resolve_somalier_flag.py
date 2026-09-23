"""
Tests for the manual flag resolution CLI.

Selection is the part worth testing. Resolution history means one key and category can match
several stored flags, so the CLI has to pick the single unresolved one or refuse, and it must never
guess. Metamist and the terminal are patched out.
"""

import pytest

from rd_qc.scripts import resolve_somalier_flag as cli

FIRST_SEEN = '2026-01-01T00:00:00+00:00'
RESOLVED_EARLIER = '2026-02-02T00:00:00+00:00'
NOW = '2026-09-22T00:00:00+00:00'
REASON = 'pedigree known wrong, family declined correction'
REVIEWER = 'ef'


def relatedness_flag(**overrides: object) -> dict:
    """A relatedness_mismatch flag as stored in SG meta, owned by CPG1 and keyed on the pair."""
    return {
        'category': 'relatedness_mismatch',
        'sequencing_group_key': 'CPG1_CPG2',
        'date': FIRST_SEEN,
        'resolved': False,
        'resolution_date': None,
        'manually_resolved': False,
        'manual_resolution_reason': None,
        'manual_resolution_by': None,
        'sg_id_1': 'CPG1',
        'sg_id_2': 'CPG2',
        'family_external_id': 'FAM1',
        'expected_relationship': 'siblings',
        'inferred_relationship': 'unrelated',
        'relatedness': 0.02,
        'ibs0': 900,
        'ibs2': 100,
    } | overrides


def held(**overrides: object) -> dict:
    """A flag already manually resolved."""
    return relatedness_flag(
        resolved=True,
        resolution_date=RESOLVED_EARLIER,
        manually_resolved=True,
        manual_resolution_reason=REASON,
        manual_resolution_by=REVIEWER,
        **overrides,
    )


def select(flags: list[dict], *, unresolve: bool = False) -> tuple[dict | None, str | None]:
    return cli.select_target(flags, 'CPG1_CPG2', 'relatedness_mismatch', 'CPG1', unresolve=unresolve)


def test_flag_key_is_order_independent():
    """A curator reads two IDs off a report row and should not have to know which sorts first."""
    assert cli.flag_key_of(['CPG2', 'CPG1']) == 'CPG1_CPG2'
    assert cli.flag_key_of(['CPG1', 'CPG2']) == 'CPG1_CPG2'
    assert cli.flag_key_of(['CPG1']) == 'CPG1'


def test_the_owning_sg_is_the_sorted_first_of_the_pair():
    """Pairwise flags are recorded against the sorted-first SG of the pair."""
    assert cli.owning_sg_id(['CPG2', 'CPG1']) == 'CPG1'
    assert cli.owning_sg_id(['CPG1']) == 'CPG1'


def test_the_one_unresolved_flag_is_selected():
    target = relatedness_flag()

    selected, problem = select([target])

    assert selected is target
    assert problem is None


def test_resolution_history_does_not_make_the_match_ambiguous():
    """The same pair flagged differently in the past is a separate, already-resolved record."""
    history = relatedness_flag(inferred_relationship='parent-child', resolved=True, resolution_date=RESOLVED_EARLIER)
    target = relatedness_flag()

    selected, _ = select([history, target])

    assert selected is target


def test_a_flag_for_another_pair_is_not_selected():
    other = relatedness_flag(sequencing_group_key='CPG1_CPG3', sg_id_2='CPG3')

    selected, problem = select([other])

    assert selected is None
    assert 'No relatedness_mismatch flag stored' in problem


def test_a_flag_of_another_category_is_not_selected():
    sex = {'category': 'sex_inference_mismatch', 'sequencing_group_key': 'CPG1', 'resolved': False}

    selected, problem = select([sex])

    assert selected is None
    assert 'No relatedness_mismatch flag stored' in problem


def test_a_legacy_flag_without_a_recorded_key_still_matches():
    """Flags written before `sequencing_group_key` existed are matched on their SG IDs instead."""
    legacy = relatedness_flag()
    del legacy['sequencing_group_key']

    selected, _ = select([legacy])

    assert selected is legacy


def test_an_already_held_flag_is_refused_with_the_reviewer_and_reason():
    selected, problem = select([held()])

    assert selected is None
    assert 'already manually resolved' in problem
    assert REVIEWER in problem
    assert '--unresolve' in problem


def test_an_auto_resolved_flag_is_refused_and_the_stored_flags_are_shown():
    """'I resolved that last month' and 'I have the wrong SG' must not read the same."""
    selected, problem = select([relatedness_flag(resolved=True, resolution_date=RESOLVED_EARLIER)])

    assert selected is None
    assert 'No unresolved relatedness_mismatch flag' in problem
    assert 'parent-child' not in problem
    assert 'siblings' in problem, 'the stored flag is printed back so the curator can see it'


def test_two_unresolved_candidates_are_refused_rather_than_guessed_at():
    """
    Reconciliation should make this unreachable, since a changed key resolves the old flag in the
    same pass. If it ever happens, picking one silently would hide a live finding.
    """
    selected, problem = select([relatedness_flag(), relatedness_flag(expected_relationship='parent-child')])

    assert selected is None
    assert 'refusing to guess' in problem


def test_unresolve_selects_the_held_flag_not_the_active_one():
    target = held()
    active = relatedness_flag(inferred_relationship='parent-child')

    selected, _ = select([target, active], unresolve=True)

    assert selected is target


def test_unresolve_is_refused_when_nothing_is_held():
    selected, problem = select([relatedness_flag()], unresolve=True)

    assert selected is None
    assert 'No manually resolved relatedness_mismatch flag' in problem


def test_with_manual_resolution_sets_the_five_fields_and_touches_nothing_else():
    original = relatedness_flag()

    updated = cli.with_manual_resolution(original, REASON, REVIEWER, NOW)

    assert updated['resolved'] is True
    assert updated['resolution_date'] == NOW
    assert updated['manually_resolved'] is True
    assert updated['manual_resolution_reason'] == REASON
    assert updated['manual_resolution_by'] == REVIEWER
    assert updated['date'] == FIRST_SEEN, 'the first-detected date is not the resolution date'
    assert updated['inferred_relationship'] == 'unrelated', 'identity is untouched'
    assert original['resolved'] is False, 'the stored flag is not mutated'


def test_without_manual_resolution_clears_every_trace():
    updated = cli.without_manual_resolution(held())

    assert updated['resolved'] is False
    assert updated['resolution_date'] is None
    assert updated['manually_resolved'] is False
    assert updated['manual_resolution_reason'] is None
    assert updated['manual_resolution_by'] is None


def test_replace_flag_swaps_one_entry_and_keeps_the_order():
    first = relatedness_flag(inferred_relationship='parent-child', resolved=True)
    target = relatedness_flag()
    replacement = cli.with_manual_resolution(target, REASON, REVIEWER, NOW)

    written = cli.replace_flag([first, target], target, replacement)

    assert written == [first, replacement]


def test_a_legacy_per_sg_flag_without_a_recorded_key_matches_on_its_owner():
    """
    Sex flags are per-SG, not pairwise, so a legacy one exercises the per-SG fallback branch of
    `sequencing_group_key` (falling back to the owning SG itself) rather than the pairwise branch
    every other test in this file goes through.
    """
    sex_flag = {
        'category': 'sex_inference_mismatch',
        'date': FIRST_SEEN,
        'resolved': False,
        'resolution_date': None,
        'manually_resolved': False,
        'manual_resolution_reason': None,
        'manual_resolution_by': None,
        'provided': 'F',
        'inferred': 'M',
        'mean_depth': 30.0,
        'x_het_ratio': 0.1,
        'x_depth_ratio': 1.0,
        'y_depth_ratio': 0.9,
        'x_sites': 500,
        'p_middling_ab': 0.02,
    }

    selected, problem = cli.select_target([sex_flag], 'CPG1', 'sex_inference_mismatch', 'CPG1', unresolve=False)

    assert selected is sex_flag
    assert problem is None


def test_replace_flag_matches_on_identity_not_equality():
    """Two history records can be equal dicts; only the selected one may be rewritten."""
    twin = relatedness_flag()
    target = relatedness_flag()
    replacement = cli.with_manual_resolution(target, REASON, REVIEWER, NOW)

    written = cli.replace_flag([twin, target], target, replacement)

    assert written[0] is twin
    assert written[1] is replacement


@pytest.fixture
def metamist(monkeypatch):
    """
    Patch the flag store, confirm the prompt, and expose what would have been written.

    `written` stays empty when nothing was written, which is what every refusal must produce.
    """
    state: dict = {'stored': [], 'written': {}, 'read_args': None}

    def fake_read(dataset: str, sg_id: str) -> list[dict]:
        state['read_args'] = (dataset, sg_id)
        return state['stored']

    def fake_write(dataset: str, sg_id: str, flags: list[dict]) -> None:
        state['written'] = {'dataset': dataset, 'sg_id': sg_id, 'flags': flags}

    monkeypatch.setattr(cli, 'read_sg_flags', fake_read)
    monkeypatch.setattr(cli, 'write_sg_flags', fake_write)
    monkeypatch.setattr(cli, 'confirmed', lambda: True)
    return state


def run(
    _metamist_state: dict,
    *,
    dataset: str = 'my-dataset',
    sg_ids: list[str] | None = None,
    category: str = 'relatedness_mismatch',
    reason: str = REASON,
    reviewer: str = REVIEWER,
    unresolve: bool = False,
    assume_yes: bool = False,
) -> int:
    """
    Invoke main with the usual arguments, overriding as needed.

    Spelled out rather than collected into a `**overrides` dict so the types survive: a dict of
    mixed values is `dict[str, object]`, which mypy cannot match against main's signature.

    `sg_ids` defaults via None because the pair is deliberately unsorted, so that every test
    exercises the sorting rather than the owning SG happening to come first.
    """
    return cli.main(
        dataset=dataset,
        sg_ids=['CPG2', 'CPG1'] if sg_ids is None else sg_ids,
        category=category,
        reason=reason,
        reviewer=reviewer,
        unresolve=unresolve,
        assume_yes=assume_yes,
    )


def test_resolving_writes_the_marked_flag_against_the_owning_sg(metamist):
    metamist['stored'] = [relatedness_flag()]

    assert run(metamist) == 0

    assert metamist['read_args'] == ('my-dataset', 'CPG1'), 'the owning SG is read, not sg_ids[0]'

    written = metamist['written']
    assert written['dataset'] == 'my-dataset'
    assert written['sg_id'] == 'CPG1', 'pairwise flags live on the sorted-first SG'
    assert len(written['flags']) == 1
    assert written['flags'][0]['manually_resolved'] is True
    assert written['flags'][0]['manual_resolution_by'] == REVIEWER
    assert written['flags'][0]['manual_resolution_reason'] == REASON
    assert written['flags'][0]['resolution_date'], 'a resolution date is stamped'


def test_the_rest_of_the_list_is_written_back_untouched(metamist):
    """The mutation replaces the whole list, so anything dropped here is lost from Metamist."""
    history = relatedness_flag(inferred_relationship='parent-child', resolved=True)
    metamist['stored'] = [history, relatedness_flag()]

    assert run(metamist) == 0

    assert metamist['written']['flags'][0] == history


def test_unresolving_clears_the_marker(metamist):
    metamist['stored'] = [held()]

    assert run(metamist, unresolve=True) == 0

    written = metamist['written']['flags'][0]
    assert written['manually_resolved'] is False
    assert written['resolved'] is False


def test_a_refused_selection_writes_nothing(metamist):
    metamist['stored'] = [held()]

    assert run(metamist) == 1
    assert metamist['written'] == {}


def test_an_sg_missing_from_the_dataset_writes_nothing(metamist, monkeypatch):
    monkeypatch.setattr(cli, 'read_sg_flags', lambda *_: None)

    assert run(metamist) == 1
    assert metamist['written'] == {}


def test_declining_the_prompt_writes_nothing(metamist, monkeypatch):
    metamist['stored'] = [relatedness_flag()]
    monkeypatch.setattr(cli, 'confirmed', lambda: False)

    assert run(metamist) == 1
    assert metamist['written'] == {}


def test_assume_yes_skips_the_prompt(metamist, monkeypatch):
    metamist['stored'] = [relatedness_flag()]

    def refuse() -> bool:
        raise AssertionError('the prompt must not be reached with --yes')

    monkeypatch.setattr(cli, 'confirmed', refuse)

    assert run(metamist, assume_yes=True) == 0


@pytest.mark.parametrize('blank', ['', '   '])
def test_an_empty_reason_or_reviewer_is_rejected_before_anything_is_read(metamist, monkeypatch, blank):
    """
    An unexplained suppression is worse than none, so this fails at the door.

    The read is patched to raise rather than merely asserting nothing was written: this CLI runs
    against production Metamist, so a malformed invocation should not cost a round trip.
    """
    metamist['stored'] = [relatedness_flag()]

    def unreachable(*_: object) -> list[dict]:
        raise AssertionError('Metamist must not be read when reason or reviewer is blank')

    monkeypatch.setattr(cli, 'read_sg_flags', unreachable)

    assert run(metamist, reason=blank) == 2
    assert run(metamist, reviewer=blank) == 2
    assert metamist['written'] == {}


def test_the_reason_and_reviewer_are_stored_stripped(metamist):
    metamist['stored'] = [relatedness_flag()]

    run(metamist, reason=f'  {REASON}  ', reviewer='  ef  ')

    written = metamist['written']['flags'][0]
    assert written['manual_resolution_reason'] == REASON
    assert written['manual_resolution_by'] == 'ef'


@pytest.mark.parametrize(
    ('typed', 'expected'),
    [
        ('y', True),
        ('Y', True),
        (' y ', True),
        ('n', False),
        ('', False),
        ('yes', True),
        ('YES', True),
        ('maybe', False),
    ],
)
def test_confirmed_accepts_only_an_explicit_affirmative(monkeypatch, typed, expected):
    """`confirmed` is the last thing standing between a curator and a production write."""
    monkeypatch.setattr('builtins.input', lambda _prompt: typed)

    assert cli.confirmed() is expected


def test_confirmed_declines_cleanly_when_stdin_is_closed(monkeypatch):
    """cron, CI, or a wrapper that forgot --yes closes stdin; that must decline, not traceback."""

    def raise_eof(_prompt: str) -> str:
        raise EOFError

    monkeypatch.setattr('builtins.input', raise_eof)

    assert cli.confirmed() is False


def test_confirmed_declines_cleanly_on_keyboard_interrupt(monkeypatch):
    """A curator hitting Ctrl-C at the prompt is a clean decline, not a traceback."""

    def raise_interrupt(_prompt: str) -> str:
        raise KeyboardInterrupt

    monkeypatch.setattr('builtins.input', raise_interrupt)

    assert cli.confirmed() is False


def test_cli_main_maps_every_argument_to_the_right_main_parameter(monkeypatch):
    """
    A typo like `unresolve=args.assume_yes` is valid Python and would ship silently.

    --unresolve is passed without --yes so the two flags disagree, which is what would catch them
    being swapped.
    """
    captured: dict = {}

    def fake_main(**kwargs: object) -> int:
        captured.update(kwargs)
        return 0

    monkeypatch.setattr(cli, 'main', fake_main)
    monkeypatch.setattr(
        'sys.argv',
        [
            'resolve_somalier_flag',
            '--dataset',
            'my-dataset',
            '--sg-ids',
            'CPG2',
            'CPG1',
            '--category',
            'relatedness_mismatch',
            '--reason',
            REASON,
            '--reviewer',
            REVIEWER,
            '--unresolve',
        ],
    )

    assert cli.cli_main() == 0
    assert captured == {
        'dataset': 'my-dataset',
        'sg_ids': ['CPG2', 'CPG1'],
        'category': 'relatedness_mismatch',
        'reason': REASON,
        'reviewer': REVIEWER,
        'unresolve': True,
        'assume_yes': False,
    }
