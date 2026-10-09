"""Command-line interface for shifted tables and Window-QROM estimates."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

from .curves import CURVES
from .tables import descriptor, materialize_table, write_manifest


def _add_public_q_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--q-x",
        help="public Q x-coordinate in hexadecimal; must be supplied with --q-y",
    )
    parser.add_argument(
        "--q-y",
        help="public Q y-coordinate in hexadecimal; must be supplied with --q-x",
    )


def _public_q(args: argparse.Namespace) -> tuple[int, int] | None:
    if (args.q_x is None) != (args.q_y is None):
        raise ValueError("--q-x and --q-y must be supplied together")
    if args.q_x is None:
        return None
    try:
        return int(args.q_x, 0), int(args.q_y, 0)
    except ValueError as error:
        raise ValueError("public Q coordinates must be hexadecimal integers") from error


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    manifest = sub.add_parser("manifest", help="write one complete (n,s) table manifest")
    manifest.add_argument("--n", type=int, choices=sorted(CURVES), required=True)
    manifest.add_argument("--s", type=int, choices=range(2, 19), required=True)
    manifest.add_argument("--output", type=Path, required=True)
    manifest.add_argument("--descriptors-only", action="store_true")
    _add_public_q_arguments(manifest)

    table = sub.add_parser("table", help="materialize one canonical raw table")
    table.add_argument("--n", type=int, choices=sorted(CURVES), required=True)
    table.add_argument("--s", type=int, choices=range(2, 19), required=True)
    table.add_argument("--stream", choices=("P", "Q"), required=True)
    table.add_argument("--j", type=int, required=True)
    table.add_argument("--basis", choices=("polynomial", "l"), required=True)
    table.add_argument("--output", type=Path, required=True)
    _add_public_q_arguments(table)

    design = sub.add_parser("design", help="synthesize one complete coherent-NCT design point")
    design.add_argument("--n", type=int, choices=sorted(CURVES), required=True)
    design.add_argument("--s", type=int, choices=range(2, 19), required=True)
    design.add_argument("--p", type=int, required=True)
    design.add_argument("--arithmetic", choices=("ac",), default="ac")
    design.add_argument("--output", type=Path, required=True)

    analytic = sub.add_parser(
        "analytic-scan",
        help="reproduce the h_u=2n module-level scan reported in the Window-QROM table",
    )
    analytic.add_argument("--results", type=Path)
    analytic.add_argument("--output-dir", type=Path, required=True)
    analytic.add_argument(
        "--n", type=int, choices=sorted(CURVES), action="append",
        help="field size to scan; repeat as needed; default scans all four",
    )
    analytic.add_argument("--s-min", type=int, default=2)
    analytic.add_argument("--s-max", type=int, default=18)

    search = sub.add_parser("search", help="run/resume the gate-exact (s,p) grid")
    search.add_argument("--n", type=int, choices=sorted(CURVES), required=True)
    search.add_argument("--arithmetic", choices=("ac",), default="ac")
    search.add_argument("--output", type=Path, required=True)
    search.add_argument("--artifact-dir", type=Path, required=True)
    search.add_argument("--s-min", type=int, default=2)
    search.add_argument("--s-max", type=int, default=18)

    coverage = sub.add_parser("coverage", help="write a fail-closed exact-grid coverage manifest")
    coverage.add_argument("--n", type=int, choices=sorted(CURVES), required=True)
    coverage.add_argument("--arithmetic", choices=("ac",), default="ac")
    coverage.add_argument("--search-csv", type=Path, required=True)
    coverage.add_argument("--output", type=Path, required=True)
    coverage.add_argument("--s-min", type=int, default=2)
    coverage.add_argument("--s-max", type=int, default=18)

    index = sub.add_parser("index", help="index completed exact-design JSON artifacts")
    index.add_argument("--n", type=int, choices=sorted(CURVES), required=True)
    index.add_argument("--arithmetic", choices=("ac",), default="ac")
    index.add_argument("--artifact-dir", type=Path, required=True)
    index.add_argument("--output", type=Path, required=True)

    tscan = sub.add_parser("toffoli-scan", help="exhaustive p=L Toffoli-count scan over s=2..18")
    tscan.add_argument("--n", type=int, choices=sorted(CURVES), required=True)
    tscan.add_argument("--arithmetic", choices=("ac",), default="ac")
    tscan.add_argument("--output", type=Path, required=True)

    verify = sub.add_parser("verify-manifest", help="verify frozen descriptors and optional table hashes")
    verify.add_argument("--manifest", type=Path, required=True)
    verify.add_argument("--table-stats", action="store_true")

    mindex = sub.add_parser("manifest-index", help="hash and index a directory of manifests")
    mindex.add_argument("--directory", type=Path, required=True)
    mindex.add_argument("--output", type=Path, required=True)
    mindex.add_argument(
        "--n",
        type=int,
        choices=sorted(CURVES),
        action="append",
        help="include only this field size; repeat to select multiple sizes",
    )

    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    if args.command == "manifest":
        manifest = write_manifest(
            args.output,
            n=args.n,
            s=args.s,
            include_table_stats=not args.descriptors_only,
            q_point=_public_q(args),
        )
        print(json.dumps({"output": str(args.output), "tables": len(manifest["tables"])}, sort_keys=True))
        return 0

    if args.command == "design":
        from .synthesis import synthesize_design, write_design

        result = synthesize_design(args.n, args.s, args.p, args.arithmetic)
        write_design(args.output, result)
        print(json.dumps(asdict(result), sort_keys=True))
        return 0

    if args.command == "analytic-scan":
        from .analytic import DEFAULT_RESULTS, FIELD_SIZES, run

        result = run(
            (args.results or DEFAULT_RESULTS).resolve(),
            args.output_dir.resolve(),
            ("upper",),
            field_sizes=args.n or FIELD_SIZES,
            s_min=args.s_min,
            s_max=args.s_max,
        )
        print(json.dumps(result, indent=2, sort_keys=True))
        return 0

    if args.command == "search":
        from .search import run_exact_grid

        run_exact_grid(
            n=args.n,
            arithmetic=args.arithmetic,
            output=args.output,
            artifact_dir=args.artifact_dir,
            s_min=args.s_min,
            s_max=args.s_max,
        )
        return 0

    if args.command == "coverage":
        from .search import write_coverage

        result = write_coverage(
            args.output,
            n=args.n,
            arithmetic=args.arithmetic,
            search_csv=args.search_csv,
            s_min=args.s_min,
            s_max=args.s_max,
        )
        print(json.dumps(result, sort_keys=True))
        return 0

    if args.command == "index":
        from .search import index_artifacts

        count = index_artifacts(
            n=args.n,
            arithmetic=args.arithmetic,
            artifact_dir=args.artifact_dir,
            output=args.output,
        )
        print(json.dumps({"exact_points": count, "output": str(args.output)}, sort_keys=True))
        return 0

    if args.command == "toffoli-scan":
        from .toffoli_search import write_toffoli_scan

        rows = write_toffoli_scan(args.output, args.n, args.arithmetic)
        best = min(rows, key=lambda row: row.total_toffoli)
        print(json.dumps(asdict(best), sort_keys=True))
        return 0

    if args.command == "verify-manifest":
        from .tables import verify_manifest

        result = verify_manifest(args.manifest, verify_table_stats=args.table_stats)
        print(json.dumps(result, sort_keys=True))
        return 0

    if args.command == "manifest-index":
        from .tables import write_manifest_index

        result = write_manifest_index(
            args.directory,
            args.output,
            field_sizes=set(args.n) if args.n else None,
        )
        print(json.dumps({"manifests": len(result["manifests"]), "output": str(args.output)}, sort_keys=True))
        return 0

    desc = descriptor(
        CURVES[args.n], args.s, args.stream, args.j, _public_q(args)
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("wb") as output:
        stats = materialize_table(
            desc,
            polynomial_output=output if args.basis == "polynomial" else None,
            l_basis_output=output if args.basis == "l" else None,
        )
    print(json.dumps({"descriptor": asdict(desc), "stats": asdict(stats), "output": str(args.output)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
