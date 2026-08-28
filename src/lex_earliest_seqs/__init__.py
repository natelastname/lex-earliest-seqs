"""Tools for defining, caching, and inspecting lexicographically earliest sequences."""

from __future__ import annotations

from .cache import (
    CacheCompatibilityError,
    cache_path_for,
    default_cache_dir,
    open_run,
)
from .core import (
    ObjectSpace,
    SequenceDefinition,
    SequenceGenerator,
    SequenceRegistry,
    SequenceRun,
    TermRecord,
)
from .incidence import (
    ColumnMode,
    IncidenceProjection,
    IncidenceRow,
    IncidenceTable,
    TextPanelMode,
    build_incidence_table,
    render_delimited,
    render_json,
    render_markdown,
    render_text,
)
from .object_space import NonNegativeIntegers, PositiveIntegers
from .zoo import register_builtins

registry = SequenceRegistry()
register_builtins(registry)

__all__ = [
    "CacheCompatibilityError",
    "ColumnMode",
    "IncidenceProjection",
    "IncidenceRow",
    "IncidenceTable",
    "NonNegativeIntegers",
    "ObjectSpace",
    "PositiveIntegers",
    "SequenceDefinition",
    "SequenceGenerator",
    "SequenceRegistry",
    "SequenceRun",
    "TermRecord",
    "TextPanelMode",
    "build_incidence_table",
    "cache_path_for",
    "default_cache_dir",
    "open_run",
    "registry",
    "render_delimited",
    "render_json",
    "render_markdown",
    "render_text",
]
