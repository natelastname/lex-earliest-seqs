"""Generic incidence chronology tables for sequence objects."""

from __future__ import annotations

import csv
import io
import json
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from enum import StrEnum
from typing import Generic, Hashable, TypeVar

from .core import TermRecord

ObjectT = TypeVar("ObjectT")
FeatureT = TypeVar("FeatureT", bound=Hashable)
CellT = TypeVar("CellT")


class ColumnMode(StrEnum):
    USED = "used"
    THROUGH_LARGEST = "through-largest"
    EXPLICIT = "explicit"


class TextPanelMode(StrEnum):
    ROWS = "rows"
    COLUMNS = "columns"


@dataclass(frozen=True)
class IncidenceProjection(Generic[ObjectT, FeatureT, CellT]):
    """Sparse coordinates used only for representation and chronology tables."""

    key: str
    title: str
    coordinates: Callable[[ObjectT], Mapping[FeatureT, CellT]]
    feature_sort_key: Callable[[FeatureT], object] = lambda value: value
    feature_label: Callable[[FeatureT], str] = str
    cell_label: Callable[[CellT], str] = str
    object_label: Callable[[ObjectT], str] = str
    through: Callable[[FeatureT], Iterable[FeatureT]] | None = None


@dataclass(frozen=True)
class IncidenceRow(Generic[ObjectT, FeatureT, CellT]):
    position: int
    subscript: int
    object_rank: int
    value: ObjectT
    coordinates: tuple[tuple[FeatureT, CellT], ...]

    def cell(self, feature: FeatureT) -> CellT | None:
        for current, value in self.coordinates:
            if current == feature:
                return value
        return None

    @property
    def support(self) -> tuple[FeatureT, ...]:
        return tuple(feature for feature, _ in self.coordinates)


@dataclass(frozen=True)
class IncidenceTable(Generic[ObjectT, FeatureT, CellT]):
    projection: IncidenceProjection[ObjectT, FeatureT, CellT]
    rows: tuple[IncidenceRow[ObjectT, FeatureT, CellT], ...]
    features: tuple[FeatureT, ...]
    omitted_features: tuple[FeatureT, ...]
    column_mode: str

    def matrix(self) -> tuple[tuple[CellT | None, ...], ...]:
        return tuple(
            tuple(row.cell(feature) for feature in self.features)
            for row in self.rows
        )


def build_incidence_table(
    records: Iterable[TermRecord[ObjectT]],
    *,
    projection: IncidenceProjection[ObjectT, FeatureT, CellT],
    column_mode: ColumnMode | str = ColumnMode.USED,
    features: Iterable[FeatureT] | None = None,
) -> IncidenceTable[ObjectT, FeatureT, CellT]:
    """Project term records into a sparse, representation-agnostic chronology table."""

    source = tuple(records)
    rows: list[IncidenceRow[ObjectT, FeatureT, CellT]] = []
    used: set[FeatureT] = set()

    for record in source:
        coordinates = projection.coordinates(record.value)
        ordered = tuple(
            sorted(
                coordinates.items(),
                key=lambda item: projection.feature_sort_key(item[0]),
            )
        )
        used.update(coordinates)
        rows.append(
            IncidenceRow(
                position=record.position,
                subscript=record.subscript,
                object_rank=record.object_rank,
                value=record.value,
                coordinates=ordered,
            )
        )

    used_features = tuple(sorted(used, key=projection.feature_sort_key))
    if features is not None:
        selected = tuple(sorted(set(features), key=projection.feature_sort_key))
        mode_name = ColumnMode.EXPLICIT.value
    else:
        mode = ColumnMode(column_mode)
        if mode is ColumnMode.EXPLICIT:
            raise ValueError("column_mode='explicit' requires features")
        if mode is ColumnMode.USED:
            selected = used_features
            mode_name = mode.value
        else:
            if projection.through is None:
                raise ValueError(
                    f"projection {projection.key!r} does not support through-largest columns"
                )
            selected = (
                tuple(projection.through(used_features[-1]))
                if used_features
                else ()
            )
            mode_name = mode.value

    selected_set = set(selected)
    omitted = tuple(feature for feature in used_features if feature not in selected_set)
    return IncidenceTable(
        projection=projection,
        rows=tuple(rows),
        features=selected,
        omitted_features=omitted,
        column_mode=mode_name,
    )


def _table_cells(table, rows, features) -> tuple[list[str], list[list[str]]]:
    projection = table.projection
    headers = [
        "n",
        "object",
        *(projection.feature_label(feature) for feature in features),
    ]
    cells: list[list[str]] = []
    for row in rows:
        values: list[str] = []
        for feature in features:
            cell = row.cell(feature)
            values.append("" if cell is None else projection.cell_label(cell))
        cells.append([str(row.subscript), projection.object_label(row.value), *values])
    return headers, cells


def _column_widths(headers, rows) -> tuple[int, ...]:
    widths = [len(header) for header in headers]
    for row in rows:
        for index, cell in enumerate(row):
            widths[index] = max(widths[index], len(cell))
    return tuple(widths)


def _render_grid(headers, rows) -> str:
    widths = _column_widths(headers, rows)

    def render_row(row) -> str:
        return " ".join(
            cell.rjust(width) for cell, width in zip(row, widths, strict=True)
        ).rstrip()

    divider = " ".join("-" * width for width in widths).rstrip()
    return "\n".join(
        [render_row(headers), divider, *(render_row(row) for row in rows)]
    )


def _rendered_width(table, rows, features) -> int:
    headers, cells = _table_cells(table, rows, features)
    return max(
        (len(line) for line in _render_grid(headers, cells).splitlines()),
        default=0,
    )


def _used_features_for_rows(table, rows):
    allowed = set(table.features)
    used = {
        feature
        for row in rows
        for feature in row.support
        if feature in allowed
    }
    return tuple(sorted(used, key=table.projection.feature_sort_key))


def _column_groups(table, rows, features, *, max_width: int):
    if not features:
        return ((),)
    if max_width <= 0:
        return (tuple(features),)

    groups = []
    current = []
    for feature in features:
        trial = tuple([*current, feature])
        if current and _rendered_width(table, rows, trial) > max_width:
            groups.append(tuple(current))
            current = [feature]
        else:
            current.append(feature)
    groups.append(tuple(current))
    return tuple(groups)


def _row_panels(table, *, max_width: int):
    if not table.rows:
        return (((), table.features),)
    if max_width <= 0:
        return ((table.rows, table.features),)

    result = []
    current = []

    def append_block(block) -> None:
        block_tuple = tuple(block)
        block_features = _used_features_for_rows(table, block_tuple)
        if _rendered_width(table, block_tuple, block_features) <= max_width:
            result.append((block_tuple, block_features))
            return
        for feature_group in _column_groups(
            table,
            block_tuple,
            block_features,
            max_width=max_width,
        ):
            result.append((block_tuple, feature_group))

    for row in table.rows:
        trial = tuple([*current, row])
        trial_features = _used_features_for_rows(table, trial)
        if current and _rendered_width(table, trial, trial_features) > max_width:
            append_block(current)
            current = [row]
        else:
            current.append(row)
    if current:
        append_block(current)
    return tuple(result)


def render_text(
    table: IncidenceTable[ObjectT, FeatureT, CellT],
    *,
    max_width: int = 120,
    panel_mode: TextPanelMode | str = TextPanelMode.ROWS,
) -> str:
    """Render a compact chronology table, preserving complete row support by default."""

    mode = TextPanelMode(panel_mode)
    use_row_panels = (
        mode is TextPanelMode.ROWS and table.column_mode == ColumnMode.USED.value
    )
    rendered: list[str] = []

    if use_row_panels:
        panels = _row_panels(table, max_width=max_width)
        for panel_index, (rows, features) in enumerate(panels, start=1):
            headers, cells = _table_cells(table, rows, features)
            body = _render_grid(headers, cells)
            if len(panels) > 1 and rows:
                body = (
                    f"rows {rows[0].subscript}–{rows[-1].subscript} "
                    f"({panel_index}/{len(panels)})\n{body}"
                )
            rendered.append(body)
    else:
        groups = _column_groups(
            table,
            table.rows,
            table.features,
            max_width=max_width,
        )
        for panel_index, features in enumerate(groups, start=1):
            headers, cells = _table_cells(table, table.rows, features)
            body = _render_grid(headers, cells)
            if len(groups) > 1:
                labels = [
                    table.projection.feature_label(feature) for feature in features
                ]
                title = (
                    f"feature columns {labels[0]}–{labels[-1]} "
                    f"({panel_index}/{len(groups)})"
                    if labels
                    else f"object columns ({panel_index}/{len(groups)})"
                )
                body = f"{title}\n{body}"
            rendered.append(body)

    if table.omitted_features:
        rendered.append(
            "omitted used features: "
            + ", ".join(
                table.projection.feature_label(feature)
                for feature in table.omitted_features
            )
        )
    return "\n\n".join(rendered) + "\n"


def _feature_chunks(features, size: int | None):
    if not features:
        return ((),)
    if size is None or size <= 0:
        return (tuple(features),)
    return tuple(
        tuple(features[offset : offset + size])
        for offset in range(0, len(features), size)
    )


def render_markdown(
    table: IncidenceTable[ObjectT, FeatureT, CellT],
    *,
    feature_columns_per_panel: int | None = 12,
) -> str:
    panels: list[str] = []
    groups = _feature_chunks(table.features, feature_columns_per_panel)
    for features in groups:
        headers, rows = _table_cells(table, table.rows, features)
        lines: list[str] = []
        if len(groups) > 1 and features:
            first = table.projection.feature_label(features[0])
            last = table.projection.feature_label(features[-1])
            lines.extend([f"### Feature columns {first}–{last}", ""])
        lines.append("| " + " | ".join(headers) + " |")
        lines.append("| " + " | ".join("---:" for _ in headers) + " |")
        lines.extend("| " + " | ".join(row) + " |" for row in rows)
        panels.append("\n".join(lines))

    if table.omitted_features:
        panels.append(
            "**Omitted used features:** "
            + ", ".join(
                table.projection.feature_label(feature)
                for feature in table.omitted_features
            )
        )
    return "\n\n".join(panels) + "\n"


def render_delimited(
    table: IncidenceTable[ObjectT, FeatureT, CellT],
    *,
    delimiter: str = ",",
) -> str:
    output = io.StringIO()
    writer = csv.writer(output, delimiter=delimiter, lineterminator="\n")
    headers, rows = _table_cells(table, table.rows, table.features)
    writer.writerow(headers)
    writer.writerows(rows)
    return output.getvalue()


def render_json(table: IncidenceTable[ObjectT, FeatureT, CellT]) -> str:
    projection = table.projection
    payload = {
        "projection": projection.key,
        "column_mode": table.column_mode,
        "features": [projection.feature_label(feature) for feature in table.features],
        "omitted_features": [
            projection.feature_label(feature) for feature in table.omitted_features
        ],
        "rows": [
            {
                "position": row.position,
                "subscript": row.subscript,
                "object_rank": row.object_rank,
                "object": projection.object_label(row.value),
                "coordinates": {
                    projection.feature_label(feature): projection.cell_label(cell)
                    for feature, cell in row.coordinates
                },
            }
            for row in table.rows
        ],
    }
    return json.dumps(payload, indent=2, sort_keys=False) + "\n"
