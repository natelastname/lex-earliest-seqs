"""Command-line interface for sequence computation and incidence tables."""

from __future__ import annotations

import argparse
from pathlib import Path

from . import registry
from .cache import open_run
from .incidence import (
    ColumnMode,
    build_incidence_table,
    render_delimited,
    render_json,
    render_markdown,
    render_text,
)


def _add_cache_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--cache-dir", type=Path, help="override the pickle cache directory"
    )
    parser.add_argument(
        "--refresh",
        action="store_true",
        help="ignore an existing pickle and regenerate",
    )
    parser.add_argument(
        "--no-cache",
        action="store_true",
        help="compute without loading or saving a pickle",
    )


def _open(args: argparse.Namespace):
    definition = registry.resolve(args.sequence)
    return open_run(
        definition,
        cache_dir=args.cache_dir,
        refresh=args.refresh,
        use_cache=not args.no_cache,
    )


def _projection(definition, name: str | None):
    projections = definition.projections or {}
    if name is None:
        if len(projections) == 1:
            return next(iter(projections.values()))
        available = ", ".join(projections) or "none"
        raise SystemExit(
            f"--projection is required; available projections: {available}"
        )
    try:
        return projections[name]
    except KeyError as exc:
        available = ", ".join(projections) or "none"
        raise SystemExit(
            f"unknown projection {name!r}; available projections: {available}"
        ) from exc


def command_list(_: argparse.Namespace) -> int:
    for definition in registry:
        oeis = f" [{definition.oeis}]" if definition.oeis else ""
        print(f"{definition.id}{oeis}\t{definition.name}")
    return 0


def command_info(args: argparse.Namespace) -> int:
    definition = registry.resolve(args.sequence)
    print(f"id: {definition.id}")
    print(f"name: {definition.name}")
    if definition.oeis:
        print(f"oeis: {definition.oeis}")
    print(f"offset: {definition.offset}")
    print(f"object space: {definition.object_space.key}")
    print(f"definition version: {definition.definition_version}")
    print(f"generator version: {definition.generator_version}")
    print("projections: " + (", ".join(definition.projections or {}) or "none"))
    if definition.description:
        print(f"description: {definition.description}")
    return 0


def command_compute(args: argparse.Namespace) -> int:
    run = _open(args)
    run.ensure(args.count)
    print(f"{run.definition.id}: cached/computed {args.count} terms")
    if run.cache_path is not None:
        print(run.cache_path)
    return 0


def command_terms(args: argparse.Namespace) -> int:
    run = _open(args)
    stop = args.start_position + args.count
    for record in run.records(args.start_position, stop):
        print(f"{record.subscript}\t{record.value}")
    return 0


def command_table(args: argparse.Namespace) -> int:
    run = _open(args)
    stop = args.start_position + args.count
    records = run.records(args.start_position, stop)
    projection = _projection(run.definition, args.projection)
    table = build_incidence_table(
        records,
        projection=projection,
        column_mode=ColumnMode(args.columns),
    )
    if args.format == "text":
        output = render_text(table, max_width=args.width)
    elif args.format == "markdown":
        output = render_markdown(table)
    elif args.format == "json":
        output = render_json(table)
    elif args.format == "csv":
        output = render_delimited(table, delimiter=",")
    else:
        output = render_delimited(table, delimiter="\t")
    print(output, end="")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="lex-earliest-seqs")
    subparsers = parser.add_subparsers(dest="command", required=True)

    list_parser = subparsers.add_parser("list", help="list built-in sequences")
    list_parser.set_defaults(func=command_list)

    info_parser = subparsers.add_parser("info", help="show sequence metadata")
    info_parser.add_argument("sequence")
    info_parser.set_defaults(func=command_info)

    compute_parser = subparsers.add_parser(
        "compute", help="compute/cache a sequence prefix"
    )
    compute_parser.add_argument("sequence")
    compute_parser.add_argument("count", type=int)
    _add_cache_options(compute_parser)
    compute_parser.set_defaults(func=command_compute)

    terms_parser = subparsers.add_parser("terms", help="print a sequence slice")
    terms_parser.add_argument("sequence")
    terms_parser.add_argument("count", type=int)
    terms_parser.add_argument("--start-position", type=int, default=0)
    _add_cache_options(terms_parser)
    terms_parser.set_defaults(func=command_terms)

    table_parser = subparsers.add_parser(
        "table", help="print an incidence chronology table"
    )
    table_parser.add_argument("sequence")
    table_parser.add_argument("count", type=int)
    table_parser.add_argument("--start-position", type=int, default=0)
    table_parser.add_argument("--projection")
    table_parser.add_argument(
        "--columns",
        choices=[ColumnMode.USED.value, ColumnMode.THROUGH_LARGEST.value],
        default=ColumnMode.USED.value,
    )
    table_parser.add_argument(
        "--format",
        choices=["text", "markdown", "json", "csv", "tsv"],
        default="text",
    )
    table_parser.add_argument("--width", type=int, default=120)
    _add_cache_options(table_parser)
    table_parser.set_defaults(func=command_table)

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    if getattr(args, "count", 0) < 0:
        parser.error("count must be nonnegative")
    if getattr(args, "start_position", 0) < 0:
        parser.error("start-position must be nonnegative")
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
