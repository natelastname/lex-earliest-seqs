"""Command-line interface for sequence computation and incidence tables."""

from __future__ import annotations

import pickle
import sys
import time
from pathlib import Path
from typing import Annotated, Literal

from cyclopts import App, Parameter, validators

from . import registry
from .cache import CacheCompatibilityError, cache_path_for, load_generator, open_run
from .incidence import (
    ColumnMode,
    build_incidence_table,
    render_delimited,
    render_json,
    render_markdown,
    render_text,
)
from .term_export import write_terms_csv, write_terms_parquet

NonNegativeInt = Annotated[
    int,
    Parameter(validator=validators.Number(gte=0)),
]
ColumnChoice = Literal["used", "through-largest"]
OutputFormat = Literal["text", "markdown", "json", "csv", "tsv"]
TermsOutputFormat = Literal["csv", "parquet"]
TermsOutputPath = Annotated[
    Path | None,
    Parameter(name=("--output", "-o")),
]

app = App(
    name="lex-earliest-seqs",
    help="Research tools for lexicographically earliest sequences.",
)


class _ProgressPrinter:
    """Throttle progress callbacks into readable stderr status lines."""

    def __init__(self, label: str, unit: Literal["terms", "bytes"]) -> None:
        self.label = label
        self.unit = unit
        self._last_time = 0.0
        self._last_value: tuple[int, int] | None = None

    @staticmethod
    def _format_bytes(value: int) -> str:
        amount = float(value)
        for suffix in ("B", "KiB", "MiB", "GiB", "TiB", "PiB"):
            if amount < 1024 or suffix == "PiB":
                if suffix == "B":
                    return f"{int(amount):,} {suffix}"
                return f"{amount:.1f} {suffix}"
            amount /= 1024
        raise AssertionError("unreachable")

    def __call__(self, current: int, total: int) -> None:
        value = (current, total)
        if value == self._last_value:
            return

        now = time.monotonic()
        done = current >= total
        if self._last_value is not None and not done and now - self._last_time < 0.5:
            return

        self._last_time = now
        self._last_value = value
        percent = 100.0 if total == 0 else min(100.0, 100.0 * current / total)
        if self.unit == "bytes":
            detail = (
                f"{self._format_bytes(current)}/{self._format_bytes(total)}"
            )
        else:
            detail = f"{current:,}/{total:,} terms"
        print(
            f"{self.label}: {detail} ({percent:5.1f}%)",
            file=sys.stderr,
            flush=True,
        )


def _open(
    sequence: str,
    *,
    cache_dir: Path | None,
    refresh: bool,
    cache: bool,
    progress: bool,
):
    definition = registry.resolve(sequence)
    load_progress = (
        _ProgressPrinter(f"load {definition.id}", "bytes") if progress else None
    )
    return open_run(
        definition,
        cache_dir=cache_dir,
        refresh=refresh,
        use_cache=cache,
        load_progress=load_progress,
    )


def _ensure(run, count: int, *, progress: bool) -> None:
    callback = (
        _ProgressPrinter(f"compute {run.definition.id}", "terms")
        if progress
        else None
    )
    run.ensure(count, progress=callback)


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


def _cached_term_count(definition, cache_dir: Path | None) -> int:
    path = cache_path_for(definition, cache_dir)
    if not path.exists():
        return 0
    return len(load_generator(definition, path).terms)


def _terms_output_format(
    output: Path | None,
    requested: TermsOutputFormat | None,
) -> TermsOutputFormat | None:
    if output is None:
        if requested is not None:
            raise SystemExit("--format requires --output/-o for the terms command")
        return None
    if requested is not None:
        return requested

    suffix = output.suffix.lower()
    if suffix == ".csv":
        return "csv"
    if suffix in {".parquet", ".pq"}:
        return "parquet"
    raise SystemExit(
        "cannot infer term output format from the filename; "
        "use --format csv or --format parquet"
    )


@app.command(name="list")
def list_sequences() -> None:
    """List built-in sequences."""

    for definition in registry:
        oeis = f" [{definition.oeis}]" if definition.oeis else ""
        print(f"{definition.id}{oeis}\t{definition.name}")


@app.command
def info(
    sequence: str,
    *,
    cache_dir: Path | None = None,
) -> None:
    """Show metadata for a sequence.

    Parameters
    ----------
    sequence
        Sequence ID, OEIS number, or registered alias.
    cache_dir
        Override the pickle cache directory when reporting the cached term count.
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
    try:
        cached_terms = _cached_term_count(definition, cache_dir)
    except (
        CacheCompatibilityError,
        EOFError,
        OSError,
        pickle.UnpicklingError,
    ) as exc:
        print(f"cached terms: unavailable ({exc})")
    else:
        print(f"cached terms: {cached_terms:,}")
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
    progress: bool = True,
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
    progress
        Print cache-loading and computation progress. Use ``--no-progress`` to
        suppress it.
    """

    run = _open(
        sequence,
        cache_dir=cache_dir,
        refresh=refresh,
        cache=cache,
        progress=progress,
    )
    _ensure(run, count, progress=progress)
    print(f"{run.definition.id}: cached/computed {count} terms")
    if run.cache_path is not None:
        print(run.cache_path)


@app.command
def terms(
    sequence: str,
    count: NonNegativeInt,
    *,
    start_position: NonNegativeInt = 0,
    output: TermsOutputPath = None,
    format: TermsOutputFormat | None = None,
    cache_dir: Path | None = None,
    refresh: bool = False,
    cache: bool = True,
    progress: bool = True,
) -> None:
    """Print or export a sequence slice.

    Parameters
    ----------
    sequence
        Sequence ID, OEIS number, or registered alias.
    count
        Number of terms to print or export.
    start_position
        Zero-based sequence position at which to start.
    output
        Write terms to this file instead of stdout. ``-o`` is an alias.
    format
        File format: ``csv`` or ``parquet``. When omitted, infer it from the
        output filename extension.
    cache_dir
        Override the pickle cache directory.
    refresh
        Ignore an existing pickle and regenerate from a fresh generator.
    cache
        Load and save the generator pickle. Use ``--no-cache`` to disable it.
    progress
        Print cache-loading and computation progress. Use ``--no-progress`` to
        suppress it.
    """

    selected_format = _terms_output_format(output, format)
    run = _open(
        sequence,
        cache_dir=cache_dir,
        refresh=refresh,
        cache=cache,
        progress=progress,
    )
    stop = start_position + count
    _ensure(run, stop, progress=progress)

    if output is None:
        for record in run.records(start_position, stop):
            print(f"{record.subscript}\t{record.value}")
    elif selected_format == "csv":
        write_terms_csv(run, start_position, stop, output)
    else:
        assert selected_format == "parquet"
        write_terms_parquet(run, start_position, stop, output)


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
    progress: bool = True,
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
    progress
        Print cache-loading and computation progress. Use ``--no-progress`` to
        suppress it.
    """

    run = _open(
        sequence,
        cache_dir=cache_dir,
        refresh=refresh,
        cache=cache,
        progress=progress,
    )
    stop = start_position + count
    _ensure(run, stop, progress=progress)
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
