"""Pickle-backed persistence for complete stateful sequence generators."""

from __future__ import annotations

import os
import pickle
import re
import tempfile
from dataclasses import dataclass, fields, is_dataclass
from pathlib import Path
from typing import BinaryIO, Generic, TypeVar

from .core import ProgressCallback, SequenceDefinition, SequenceGenerator, SequenceRun

ObjectT = TypeVar("ObjectT")
CACHE_FORMAT_VERSION = 2
GeneratorSchema = tuple[str, str, str, tuple[tuple[str, str], ...]]


class CacheCompatibilityError(RuntimeError):
    """Raised when a pickle cache belongs to incompatible code or metadata."""


@dataclass
class CachedGenerator(Generic[ObjectT]):
    cache_format: int
    definition_id: str
    definition_version: int
    generator_version: int
    generator_schema: GeneratorSchema
    generator: SequenceGenerator[ObjectT]


class _ProgressReader:
    """Binary-file proxy that reports bytes consumed by ``pickle.load``."""

    def __init__(
        self,
        handle: BinaryIO,
        total_bytes: int,
        progress: ProgressCallback,
    ) -> None:
        self._handle = handle
        self._total_bytes = total_bytes
        self._progress = progress

    def _report(self) -> None:
        self._progress(min(self._handle.tell(), self._total_bytes), self._total_bytes)

    def read(self, size: int = -1) -> bytes:
        data = self._handle.read(size)
        self._report()
        return data

    def readline(self, size: int = -1) -> bytes:
        data = self._handle.readline(size)
        self._report()
        return data

    def readinto(self, buffer) -> int | None:
        count = self._handle.readinto(buffer)
        self._report()
        return count


def default_cache_dir() -> Path:
    root = os.environ.get("XDG_CACHE_HOME")
    if root:
        return Path(root) / "lex-earliest-seqs"
    return Path.home() / ".cache" / "lex-earliest-seqs"


def cache_filename(definition: SequenceDefinition[object]) -> str:
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "-", definition.id).strip("-.")
    if not safe:
        raise ValueError(f"sequence id cannot form a cache filename: {definition.id!r}")
    return f"{safe}.pkl"


def cache_path_for(
    definition: SequenceDefinition[object],
    cache_dir: str | os.PathLike[str] | None = None,
) -> Path:
    directory = Path(cache_dir) if cache_dir is not None else default_cache_dir()
    return directory / cache_filename(definition)


def _annotation_name(annotation: object) -> str:
    if isinstance(annotation, str):
        return annotation
    module = getattr(annotation, "__module__", None)
    qualname = getattr(annotation, "__qualname__", None)
    if module is not None and qualname is not None:
        return f"{module}.{qualname}"
    return repr(annotation)


def _generator_schema(generator: SequenceGenerator[object]) -> GeneratorSchema:
    """Return a stable fingerprint of the generator's persisted state layout."""

    generator_type = type(generator)
    if is_dataclass(generator):
        kind = "dataclass"
        state_fields = tuple(
            (field.name, _annotation_name(field.type)) for field in fields(generator)
        )
    else:
        kind = "annotations"
        annotations: dict[str, object] = {}
        for base in reversed(generator_type.__mro__):
            annotations.update(getattr(base, "__annotations__", {}))
        state_fields = tuple(
            sorted(
                (name, _annotation_name(annotation))
                for name, annotation in annotations.items()
            )
        )
    return (
        generator_type.__module__,
        generator_type.__qualname__,
        kind,
        state_fields,
    )


def _fresh_generator(
    definition: SequenceDefinition[ObjectT],
) -> SequenceGenerator[ObjectT]:
    generator = definition.generator_factory()
    if not isinstance(generator, SequenceGenerator):
        raise TypeError("generator_factory did not return a SequenceGenerator")
    return generator


def _validate_cached(
    definition: SequenceDefinition[ObjectT],
    cached: object,
    expected_schema: GeneratorSchema,
) -> CachedGenerator[ObjectT]:
    if not isinstance(cached, CachedGenerator):
        raise CacheCompatibilityError("cache does not contain a CachedGenerator")
    expected = (
        CACHE_FORMAT_VERSION,
        definition.id,
        definition.definition_version,
        definition.generator_version,
    )
    actual = (
        cached.cache_format,
        cached.definition_id,
        cached.definition_version,
        cached.generator_version,
    )
    if actual != expected:
        raise CacheCompatibilityError(
            "incompatible sequence cache: "
            f"expected {expected!r}, found {actual!r}"
        )
    actual_schema = getattr(cached, "generator_schema", None)
    if actual_schema != expected_schema:
        raise CacheCompatibilityError(
            "incompatible generator cache schema: "
            f"expected {expected_schema!r}, found {actual_schema!r}"
        )
    if not isinstance(cached.generator, SequenceGenerator):
        raise CacheCompatibilityError("cached generator does not satisfy SequenceGenerator")
    return cached


def load_generator(
    definition: SequenceDefinition[ObjectT],
    path: str | os.PathLike[str],
    *,
    progress: ProgressCallback | None = None,
) -> SequenceGenerator[ObjectT]:
    cache_path = Path(path)
    expected_schema = _generator_schema(_fresh_generator(definition))
    total_bytes = cache_path.stat().st_size
    try:
        with cache_path.open("rb") as handle:
            try:
                if progress is None:
                    cached = pickle.load(handle)
                else:
                    progress(0, total_bytes)
                    cached = pickle.load(_ProgressReader(handle, total_bytes, progress))
                    progress(total_bytes, total_bytes)
            except (AttributeError, ImportError) as exc:
                raise CacheCompatibilityError(
                    "cache cannot be loaded with the current generator schema"
                ) from exc
        return _validate_cached(definition, cached, expected_schema).generator
    except CacheCompatibilityError:
        cache_path.unlink(missing_ok=True)
        raise


def save_generator(
    definition: SequenceDefinition[ObjectT],
    generator: SequenceGenerator[ObjectT],
    path: str | os.PathLike[str],
) -> Path:
    cache_path = Path(path)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cached = CachedGenerator(
        cache_format=CACHE_FORMAT_VERSION,
        definition_id=definition.id,
        definition_version=definition.definition_version,
        generator_version=definition.generator_version,
        generator_schema=_generator_schema(generator),
        generator=generator,
    )

    fd, temporary_name = tempfile.mkstemp(
        dir=cache_path.parent,
        prefix=f".{cache_path.name}.",
        suffix=".tmp",
    )
    temporary_path = Path(temporary_name)
    try:
        with os.fdopen(fd, "wb") as handle:
            pickle.dump(cached, handle, protocol=pickle.HIGHEST_PROTOCOL)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_path, cache_path)
    except BaseException:
        temporary_path.unlink(missing_ok=True)
        raise
    return cache_path


def open_run(
    definition: SequenceDefinition[ObjectT],
    *,
    cache_dir: str | os.PathLike[str] | None = None,
    cache_path: str | os.PathLike[str] | None = None,
    refresh: bool = False,
    use_cache: bool = True,
    load_progress: ProgressCallback | None = None,
) -> SequenceRun[ObjectT]:
    if cache_dir is not None and cache_path is not None:
        raise ValueError("provide cache_dir or cache_path, not both")

    path: Path | None
    if use_cache:
        path = (
            Path(cache_path)
            if cache_path is not None
            else cache_path_for(definition, cache_dir)
        )
    else:
        path = None

    if path is not None and path.exists() and not refresh:
        try:
            generator = load_generator(definition, path, progress=load_progress)
        except CacheCompatibilityError:
            generator = _fresh_generator(definition)
    else:
        generator = _fresh_generator(definition)

    return SequenceRun(definition, generator, cache_path=path)
