"""
Render the Somalier flags report from mock fixtures, with no Metamist access at all.

Deterministic, and it dodges the transient Metamist GraphQL stalls, so this is the way to iterate
on the page before trying a live run. Writes two files into the output directory (the system temp
directory by default):

    somalier_flags_report.html            the full fixture set
    somalier_flags_report_all_clear.html  the same dataset with no flags

Run with:  uv run python testing_scripts/render_mock_somalier_report.py
"""

import sys
from argparse import ArgumentParser
from pathlib import Path
from tempfile import gettempdir

REPO_ROOT = Path(__file__).absolute().parent.parent
# The fixtures live under test/ so pytest can import them with no path setup. This dev script gets
# at them the same way pytest does, by putting test/ on the path. src/ goes on too, because a local
# .venv often holds a stale non-editable rd_qc (see the pythonpath note in pyproject.toml).
sys.path.insert(0, str(REPO_ROOT / 'test'))
sys.path.insert(0, str(REPO_ROOT / 'src'))

from fixtures.somalier_flags import (  # noqa: E402
    MOCK_ALL_CLEAR_SEQUENCING_GROUPS,
    MOCK_SEQUENCING_GROUPS,
    MOCK_SG_INFOS,
)

from rd_qc.scripts.somalier_flags_report import (  # noqa: E402
    FamilyGroup,
    SGInfo,
    collect_somalier_flags,
    group_by_family,
    render_report,
    split_active_resolved,
    summarise_flags,
)

DATASET = 'mock-relatedness-test'
GENERATED_AT = '2026-09-14T10:00:00+00:00'


def build_html(
    sequencing_groups: list[dict],
    infos: dict[str, SGInfo],
) -> tuple[str, dict, list[FamilyGroup], list[FamilyGroup]]:
    """Run the fixture set through the whole report pipeline and render it."""
    flagged = [sf for sf in collect_somalier_flags(sequencing_groups) if sf.flags]
    groups = group_by_family(flagged, infos)
    active_groups, resolved_groups = split_active_resolved(groups, infos)
    summary = summarise_flags(
        flagged,
        total_sgs=len(sequencing_groups),
        families_affected=len(active_groups),
        infos=infos,
    )
    html = render_report(
        DATASET,
        active_groups,
        resolved_groups,
        summary=summary,
        generated_at=GENERATED_AT,
    )
    return html, summary, active_groups, resolved_groups


def describe(
    label: str,
    summary: dict,
    active_groups: list[FamilyGroup],
    resolved_groups: list[FamilyGroup],
) -> None:
    """Print what got rendered, so the page can be sanity checked without opening a browser."""
    print(f'\n{label}')
    print(f'  summary: {summary}')
    print(f'  active groups ({len(active_groups)}):')
    for group in active_groups:
        print(f'    {group.label:<26} {group.count_summary}')
        for flag in group.flags:
            cross = f'  [cross-family -> {flag.cross_family}]' if flag.cross_family else ''
            print(f'        {flag.category_key:<9} {flag.subject:<24} {flag.result}{cross}')
    print(f'  resolved groups ({len(resolved_groups)}):')
    for group in resolved_groups:
        print(f'    {group.label:<26} {group.count_summary}')
        for flag in group.flags:
            print(f'        {flag.category_key:<9} {flag.subject:<24} {flag.result}')


def render_to(output_dir: Path, name: str, sequencing_groups: list[dict]) -> Path:
    html, summary, active_groups, resolved_groups = build_html(sequencing_groups, MOCK_SG_INFOS)
    out = output_dir / name
    out.write_text(html)
    describe(f'Wrote {out} ({len(html):,} bytes)', summary, active_groups, resolved_groups)
    return out


def main(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    render_to(output_dir, 'somalier_flags_report.html', MOCK_SEQUENCING_GROUPS)
    render_to(output_dir, 'somalier_flags_report_all_clear.html', MOCK_ALL_CLEAR_SEQUENCING_GROUPS)


if __name__ == '__main__':
    parser = ArgumentParser(description='Render the Somalier flags report from mock fixtures.')
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path(gettempdir()),
        help='directory to write the rendered HTML into (default: the system temp directory)',
    )
    main(parser.parse_args().output_dir)
