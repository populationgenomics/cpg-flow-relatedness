"""
Record or clear a manual resolution on a single Somalier flag in Metamist.

A finding a curator has reviewed and accepted (a pedigree known to be wrong, a mismatch already
investigated and explained) reappears in the report every run, because reconciliation reopens any
flag it measures again. Marking it manually resolved holds it resolved across runs and keeps it out
of the report, without stopping the checks from measuring it.

Run locally, against your own Metamist credentials:

    resolve_somalier_flag --dataset my-dataset --sg-ids CPG001 CPG002 \
        --category relatedness_mismatch --reason "pedigree known wrong" --reviewer ef

Take the SG IDs off the report row in either order. Add --unresolve to reopen a flag.
"""

import json
from argparse import ArgumentParser
from datetime import UTC, datetime

from loguru import logger

from rd_qc.flag_store import read_sg_flags, sequencing_group_key, write_sg_flags
from rd_qc.utils import sg_ids_tag

CATEGORIES = ('sex_inference_mismatch', 'self_relatedness_mismatch', 'relatedness_mismatch')

EXIT_OK = 0
EXIT_REFUSED = 1
EXIT_BAD_ARGS = 2


def flag_key_of(sg_ids: list[str]) -> str:
    """The `sequencing_group_key` for these SG IDs, given in either order."""
    return sg_ids_tag(sg_ids)


def owning_sg_id(sg_key: str) -> str:
    """
    The sequencing group whose meta holds a flag with this key.

    Pairwise flags are recorded against the sorted-first SG of the pair and the key is the sorted
    join, so the owner is the first element.
    """
    return sg_key.split('_', maxsplit=1)[0]


def describe(flags: list[dict]) -> str:
    """One indented JSON line per flag, for printing candidates back to the curator."""
    return '\n'.join(f'  {json.dumps(flag, sort_keys=True)}' for flag in flags)


def matching_flags(flags: list[dict], sg_key: str, category: str, owner: str) -> list[dict]:
    """
    Every stored flag for this key and category, resolved or not.

    Falls back to deriving the key for flags recorded before `sequencing_group_key` existed, the
    same way the report does.
    """
    return [
        flag
        for flag in flags
        if flag.get('category') == category
        and (flag.get('sequencing_group_key') or sequencing_group_key(flag, owner)) == sg_key
    ]


def select_target(
    flags: list[dict],
    sg_key: str,
    category: str,
    owner: str,
    *,
    unresolve: bool,
) -> tuple[dict | None, str]:
    """
    The single flag to act on, or `None` and the reason why not.

    Resolving acts on an unresolved flag, unresolving on a manually resolved one. Anything other
    than exactly one candidate is refused rather than guessed at: resolution history means a key
    and category can match several stored records, and picking one silently could hide a live
    finding. Refusals print the stored flags, so 'I resolved that last month' and 'I have the wrong
    sequencing group' do not read the same.
    """
    stored = matching_flags(flags, sg_key, category, owner)
    if not stored:
        return None, f'No {category} flag stored for {sg_key}. Check the SG IDs and the category.'

    if unresolve:
        candidates = [flag for flag in stored if flag.get('manually_resolved')]
        if not candidates:
            return None, f'No manually resolved {category} flag for {sg_key}. Stored flags:\n{describe(stored)}'
    else:
        candidates = [flag for flag in stored if not flag.get('resolved', False)]
        if not candidates:
            already_held = [flag for flag in stored if flag.get('manually_resolved')]
            if already_held:
                flag = already_held[0]
                return None, (
                    f'{category} for {sg_key} is already manually resolved, on '
                    f'{flag.get("resolution_date")} by {flag.get("manual_resolution_by")}: '
                    f'{flag.get("manual_resolution_reason")}. Use --unresolve to reopen it.'
                )
            return None, f'No unresolved {category} flag for {sg_key}. Stored flags:\n{describe(stored)}'

    if len(candidates) > 1:
        return None, (
            f'{len(candidates)} candidate {category} flags for {sg_key}; refusing to guess.\n{describe(candidates)}'
        )
    return candidates[0], ''


def with_manual_resolution(flag: dict, reason: str, reviewer: str, now: str) -> dict:
    """
    A copy of `flag`, resolved by hand.

    `resolved` and `resolution_date` are set alongside the manual fields so a reader that knows
    nothing about manual resolution still sees a resolved flag. Identity fields and the
    first-detected `date` are untouched.
    """
    return flag | {
        'resolved': True,
        'resolution_date': now,
        'manually_resolved': True,
        'manual_resolution_reason': reason,
        'manual_resolution_by': reviewer,
    }


def without_manual_resolution(flag: dict) -> dict:
    """A copy of `flag`, active again, with every trace of the manual resolution cleared."""
    return flag | {
        'resolved': False,
        'resolution_date': None,
        'manually_resolved': False,
        'manual_resolution_reason': None,
        'manual_resolution_by': None,
    }


def replace_flag(flags: list[dict], target: dict, replacement: dict) -> list[dict]:
    """
    A new list with `target` swapped for `replacement`, matched on identity rather than equality.

    Two history records for one pair can be equal dicts, and only the selected one may be rewritten.
    """
    return [replacement if flag is target else flag for flag in flags]
