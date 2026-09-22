"""
Tests for the one module that reads and writes SG meta's `somalier_flags` list.

The mutation overwrites the whole list, so both the reconciler and the resolve CLI go through here.
Metamist is patched out: what matters is the shape of what would be sent.
"""

from rd_qc import flag_store

FLAG = {'category': 'relatedness_mismatch', 'sg_id_1': 'CPG1', 'sg_id_2': 'CPG2', 'resolved': False}


def fake_dataset(monkeypatch, sequencing_groups: list[dict]) -> dict:
    """
    Patch `query` to answer the dataset read, and capture any mutation variables.

    This ignores which GraphQL document it was handed, so it answers every query with the same
    read-shaped payload. Fine while there is one read and one mutation; a second read of a
    different shape would get the wrong answer here rather than a loud failure.
    """
    captured: dict = {}

    def fake_query(_query, variables=None) -> dict:
        captured.update(variables or {})
        return {'project': {'sequencingGroups': sequencing_groups}}

    monkeypatch.setattr(flag_store, 'query', fake_query)
    return captured


def test_read_sg_flags_returns_the_stored_list(monkeypatch):
    fake_dataset(monkeypatch, [{'id': 'CPG1', 'meta': {'somalier_flags': [FLAG]}}])

    assert flag_store.read_sg_flags('my-dataset', 'CPG1') == [FLAG]


def test_read_sg_flags_returns_empty_for_a_sequencing_group_with_no_flags(monkeypatch):
    """
    An SG that exists but has never been flagged is not the same as a missing SG.

    Its meta is populated with other things, so this is the common real shape: the key is absent
    rather than the meta being empty.
    """
    fake_dataset(monkeypatch, [{'id': 'CPG1', 'meta': {'sequencing_type': 'genome'}}])

    assert flag_store.read_sg_flags('my-dataset', 'CPG1') == []


def test_read_sg_flags_returns_empty_for_a_sequencing_group_with_null_meta(monkeypatch):
    """Metamist returns meta as null rather than {} for some SGs, which must not raise."""
    fake_dataset(monkeypatch, [{'id': 'CPG1', 'meta': None}])

    assert flag_store.read_sg_flags('my-dataset', 'CPG1') == []


def test_read_sg_flags_returns_none_for_a_sequencing_group_not_in_the_dataset(monkeypatch):
    """Distinguishable from the empty case so a caller can say 'wrong dataset' rather than 'no flags'."""
    fake_dataset(monkeypatch, [{'id': 'CPG1', 'meta': {}}])

    assert flag_store.read_sg_flags('my-dataset', 'CPG9') is None


def test_write_sg_flags_sends_the_whole_list_under_the_meta_key(monkeypatch):
    captured = fake_dataset(monkeypatch, [])

    flag_store.write_sg_flags('my-dataset', 'CPG1', [FLAG])

    assert captured == {'dataset': 'my-dataset', 'sgId': 'CPG1', 'sgMeta': {'somalier_flags': [FLAG]}}


def test_sequencing_group_key_joins_a_pair_in_sorted_order():
    assert flag_store.sequencing_group_key({'sg_id_1': 'CPG2', 'sg_id_2': 'CPG1'}, 'CPG2') == 'CPG1_CPG2'


def test_sequencing_group_key_falls_back_to_the_owning_sg_for_a_per_sg_flag():
    assert flag_store.sequencing_group_key({'provided': 'M', 'inferred': 'F'}, 'CPG1') == 'CPG1'
