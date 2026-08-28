"""Command-line interface for sequence computation and incidence tables."""

from __future__ import annotations

from pathlib import Path
from typing import Annotated, Literal

from cyclopts import App, Parameter, validators

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

NonNegativeInt = Annotated[
    int,
    Parameter(validator=validators.Number(gte=0)),
]
ColumnChoice = Literal["used", "through-largest"]
OutputFormat = Literal["text", "markdown", "json", "csv", "tsv"]

app = App(
    name="lex-earliest-seqs",
    help="Research tools for lexicographically earliest sequences.",
)


def _open(
    sequence: str,
    *,
    cache_dir: Path | None,
    refresh: bool,
    cache: bool,
):
    definition = registry.resolve(sequence)
    return open_run(
        definition,
        cache_dir=cache_dir,
        refresh=refresh,
        use_cache=cache,
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


@app.command(name="list")
def list_sequences() -> None:
    """List built-in sequences."""

    for definition in registry:
        oeis = f" [{definition.oeis}]" if definition.oeis else ""
        print(f"{definition.id}{oeis}\t{definition.name}")


@app.command
def info(sequence: str) -> None:
    """Show metadata for a sequence.

    Parameters
    ----------
    sequence
        Sequence ID, OEIS number, or registered alias.
    """

    definition = registry.resolve(sequence)
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


@app.command
def compute(
    sequence: str,
    count: NonNegativeInt,
    *,
    cache_dir: Path | None = None,
    refresh: bool = False,
    cache: bool = True,
) -> None:
    """Compute and optionally cache a sequence prefix.

    Parameters
    ----------
    sequence
        Sequence ID, OEIS number, or registered alias.
    count
        Number of terms to ensure are available.
    cache_dir
        Override the pickle cache directory.
    refresh
        Ignore an existing pickle and regenerate from a fresh generator.
    cache
        Load and save the generator pickle. Use ``--no-cache`` to disable it.
    """

    run = _open(
        sequence,
        cache_dir=cache_dir,
        refresh=refresh,
        cache=cache,
    )
    run.ensure(count)
    print(f"{run.definition.id}: cached/computed {count} terms")
    if run.cache_path is not None:
        print(run.cache_path)


@app.command
def terms(
    sequence: str,
    count: NonNegativeInt,
    *,
    start_position: NonNegativeInt = 0,
    cache_dir: Path | None = None,
    refresh: bool = False,
    cache: bool = True,
) -> None:
    """Print a sequence slice.

    Parameters
    ----------
    sequence
        Sequence ID, OEIS number, or registered alias.
    count
        Number of terms to print.
    start_position
        Zero-based sequence position at which to start.
    cache_dir
        Override the pickle cache directory.
    refresh
        Ignore an existing pickle and regenerate from a fresh generator.
    cache
        Load and save the generator pickle. Use ``--no-cache`` to disable it.
    """

    run = _open(
        sequence,
        cache_dir=cache_dir,
        refresh=refresh,
        cache=cache,
    )
    stop = start_position + count
    for record in run.records(start_position, stop):
        print(f"{record.subscript}\t{record.value}")


@app.command
def table(
    sequence: str,
    count: NonNegativeInt,
    *,
    start_position: NonNegativeInt = 0,
    projection: str | None = None,
    columns: ColumnChoice = "used",
    format: OutputFormat = "text",
    width: NonNegativeInt = 120,
    cache_dir: Path | None = None,
    refresh: bool = False,
    cache: bool = True,
) -> None:
    """Print an incidence chronology table.

    Parameters
    ----------
    sequence
        Sequence ID, OEIS number, or registered alias.
    count
        Number of sequence terms to include.
    start_position
        Zero-based sequence position at which to start.
    projection
        Named incidence projection. Omit when the sequence has exactly one.
    columns
        Feature-column policy: ``used`` or ``through-largest``.
    format
        Output format: text, markdown, json, csv, or tsv.
    width
        Maximum text-table width before panel splitting.
    cache_dir
        Override the pickle cache directory.
    refresh
        Ignore an existing pickle and regenerate from a fresh generator.
    cache
        Load and save the generator pickle. Use ``--no-cache`` to disable it.
    """

    run = _open(
        sequence,
        cache_dir=cache_dir,
        refresh=refresh,
        cache=cache,
    )
    stop = start_position + count
    records = run.records(start_position, stop)
    selected_projection = _projection(run.definition, projection)
    incidence_table = build_incidence_table(
        records,
        projection=selected_projection,
        column_mode=ColumnMode(columns),
    )

    if format == "text":
        output = render_text(incidence_table, max_width=width)
    elif format == "markdown":
        output = render_markdown(incidence_table)
    elif format == "json":
        output = render_json(incidence_table)
    elif format == "csv":
        output = render_delimited(incidence_table, delimiter=",")
    else:
        output = render_delimited(incidence_table, delimiter="\t")
    print(output, end="")


def main() -> None:
    """Run the Cyclopts application."""

    app()


if __name__ == "__main__":
    main()
