"""CSV and Parquet export for sequence terms."""

from __future__ import annotations

import csv
import os
import tempfile
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq

from .core import SequenceRun

DEFAULT_EXPORT_BATCH_SIZE = 131_072


def _temporary_output_path(output: Path) -> Path:
    output.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(
        dir=output.parent,
        prefix=f".{output.name}.",
        suffix=".tmp",
    )
    os.close(fd)
    return Path(name)


def _validate_range(start: int, stop: int, batch_size: int) -> None:
    if start < 0:
        raise ValueError("start must be nonnegative")
    if stop < start:
        raise ValueError("stop must be at least start")
    if batch_size < 1:
        raise ValueError("batch_size must be positive")


def write_terms_csv(
    run: SequenceRun[Any],
    start: int,
    stop: int,
    output: str | Path,
    *,
    batch_size: int = DEFAULT_EXPORT_BATCH_SIZE,
) -> Path:
    """Write ``subscript,value`` rows to CSV without an extra full-prefix copy."""

    _validate_range(start, stop, batch_size)
    run.ensure(stop)
    path = Path(output)
    temporary = _temporary_output_path(path)
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(("subscript", "value"))
            for batch_start in range(start, stop, batch_size):
                batch_stop = min(stop, batch_start + batch_size)
                values = run.terms[batch_start:batch_stop]
                writer.writerows(
                    (position + run.definition.offset, value)
                    for position, value in zip(
                        range(batch_start, batch_stop), values, strict=True
                    )
                )
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    return path


def _term_table(run: SequenceRun[Any], start: int, stop: int) -> pa.Table:
    values = run.terms[start:stop]
    try:
        value_array = pa.array(values)
    except (pa.ArrowException, OverflowError, TypeError) as exc:
        raise TypeError(
            "Parquet export requires term values that PyArrow can represent "
            "as a scalar column."
        ) from exc
    return pa.table(
        {
            "subscript": pa.array(
                range(
                    start + run.definition.offset,
                    stop + run.definition.offset,
                ),
                type=pa.int64(),
            ),
            "value": value_array,
        }
    )


def write_terms_parquet(
    run: SequenceRun[Any],
    start: int,
    stop: int,
    output: str | Path,
    *,
    batch_size: int = DEFAULT_EXPORT_BATCH_SIZE,
) -> Path:
    """Write ``subscript,value`` rows to compressed Parquet in bounded batches."""

    _validate_range(start, stop, batch_size)
    run.ensure(stop)
    path = Path(output)
    temporary = _temporary_output_path(path)
    writer: pq.ParquetWriter | None = None
    try:
        if start == stop:
            empty = pa.table(
                {
                    "subscript": pa.array([], type=pa.int64()),
                    "value": pa.array([], type=pa.null()),
                }
            )
            pq.write_table(empty, temporary, compression="zstd")
        else:
            for batch_start in range(start, stop, batch_size):
                batch_stop = min(stop, batch_start + batch_size)
                table = _term_table(run, batch_start, batch_stop)
                if writer is None:
                    writer = pq.ParquetWriter(
                        temporary,
                        table.schema,
                        compression="zstd",
                    )
                elif table.schema != writer.schema:
                    try:
                        table = table.cast(writer.schema)
                    except (pa.ArrowException, ValueError) as exc:
                        raise TypeError(
                            "term values changed Arrow type between export batches"
                        ) from exc
                writer.write_table(table)
            assert writer is not None
            writer.close()
            writer = None
        os.replace(temporary, path)
    except BaseException:
        if writer is not None:
            writer.close()
        temporary.unlink(missing_ok=True)
        raise
    return path
