"""Temporary migration helper for the legacy enots-wolley-2 EW cache.

This module is intentionally narrow and disposable. It converts the historical
list-only A336957 pickle into the native stateful-generator cache without
replaying any greedy sequence decisions.
"""

from __future__ import annotations

import gc
import pickle
import sys
from pathlib import Path
from typing import BinaryIO

from . import registry
from .cache import cache_path_for, load_generator, save_generator
from .zoo.enots_wolley import EnotsWolleyGenerator

LEGACY_FORMAT_VERSION = 1
LEGACY_SEQUENCE_ID = "oeis-a336957-original-v1"
DEFAULT_LEGACY_CACHE_PATH = (
    Path.home() / ".cache" / "enots-wolley-2" / "terms-v1.pkl"
)
EXPECTED_PREFIX = [1, 2, 6, 15, 35, 14, 12, 33, 55, 10, 18, 21]


class LegacyEWCacheError(RuntimeError):
    """Raised when the source pickle is not the expected historical EW cache."""


class _ProgressReader:
    """Minimal pickle reader that reports roughly twenty byte checkpoints."""

    def __init__(self, handle: BinaryIO, total_bytes: int, label: str) -> None:
        self._handle = handle
        self._total_bytes = total_bytes
        self._label = label
        self._step = max(1, total_bytes // 20)
        self._next_report = 0

    def _report(self) -> None:
        consumed = min(self._handle.tell(), self._total_bytes)
        if consumed < self._next_report and consumed != self._total_bytes:
            return
        percentage = 100.0 if self._total_bytes == 0 else 100.0 * consumed / self._total_bytes
        print(
            f"{self._label}: {consumed:,}/{self._total_bytes:,} bytes "
            f"({percentage:5.1f}%)",
            file=sys.stderr,
            flush=True,
        )
        while self._next_report <= consumed:
            self._next_report += self._step

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


def _load_legacy_terms(source: Path, *, progress: bool) -> list[int]:
    if not source.exists():
        raise FileNotFoundError(source)

    total_bytes = source.stat().st_size
    if progress:
        print(f"loading legacy cache: {source} ({total_bytes:,} bytes)", file=sys.stderr)

    with source.open("rb") as handle:
        payload = pickle.load(
            _ProgressReader(handle, total_bytes, "legacy load") if progress else handle
        )

    if not isinstance(payload, dict):
        raise LegacyEWCacheError("legacy cache payload is not a dict")
    if payload.get("format_version") != LEGACY_FORMAT_VERSION:
        raise LegacyEWCacheError(
            f"unexpected legacy format_version: {payload.get('format_version')!r}"
        )
    if payload.get("sequence_id") != LEGACY_SEQUENCE_ID:
        raise LegacyEWCacheError(
            f"unexpected legacy sequence_id: {payload.get('sequence_id')!r}"
        )

    terms = payload.get("terms")
    if not isinstance(terms, list):
        raise LegacyEWCacheError("legacy cache terms are not a list")
    if payload.get("term_count") != len(terms):
        raise LegacyEWCacheError(
            "legacy term_count does not match the stored term list length"
        )
    if terms[: len(EXPECTED_PREFIX)] != EXPECTED_PREFIX[: len(terms)]:
        raise LegacyEWCacheError("legacy cache does not have the expected A336957 prefix")
    if any(type(term) is not int or term < 1 for term in terms):
        raise LegacyEWCacheError("legacy cache contains a non-positive-integer term")

    return terms


def _next_power_of_two_at_least(value: int, *, minimum: int = 1_024) -> int:
    limit = minimum
    while limit < value:
        limit *= 2
    return limit


def _generator_from_legacy_terms(
    terms: list[int],
    *,
    progress: bool,
) -> EnotsWolleyGenerator:
    if progress:
        print(f"reconstructing used set from {len(terms):,} terms", file=sys.stderr)
    used = set(terms)
    if len(used) != len(terms):
        raise LegacyEWCacheError("legacy cache contains duplicate terms")

    smallest_unused = 1
    while smallest_unused in used:
        smallest_unused += 1

    # Only the predecessor and the next scan point need to fit in the radical
    # table when computation resumes. The table itself is deliberately absent
    # from the migrated pickle and will be generated lazily on first extension.
    needed = smallest_unused
    if terms:
        needed = max(needed, terms[-1])
    if len(terms) >= 2:
        needed = max(needed, terms[-2])
    limit = _next_power_of_two_at_least(needed)

    if progress:
        print(
            f"reconstructed state: smallest_unused={smallest_unused:,}, "
            f"radical capacity={limit:,}",
            file=sys.stderr,
        )

    return EnotsWolleyGenerator(
        terms=terms,
        used=used,
        smallest_unused=smallest_unused,
        limit=limit,
        radicals=None,
    )


def migrate_legacy_ew_cache(
    source: Path = DEFAULT_LEGACY_CACHE_PATH,
    destination: Path | None = None,
    *,
    force: bool = False,
    verify: bool = True,
    progress: bool = True,
) -> Path:
    """Convert the historical A336957 prefix cache into the native pickle.

    No sequence terms are generated. Continuation state is reconstructed solely
    from the trusted legacy prefix.
    """

    definition = registry.resolve("A336957")
    target = destination or cache_path_for(definition)
    source = source.expanduser()
    target = target.expanduser()

    if target.exists() and not force:
        raise FileExistsError(
            f"destination already exists: {target}; pass --force to replace it"
        )

    terms = _load_legacy_terms(source, progress=progress)
    generator = _generator_from_legacy_terms(terms, progress=progress)
    term_count = len(generator.terms)
    last_term = generator.terms[-1] if generator.terms else None
    smallest_unused = generator.smallest_unused

    if progress:
        print(f"writing native cache: {target}", file=sys.stderr)
    save_generator(definition, generator, target)

    # Drop the original 5M-term structures before verification so the second
    # unpickle does not coexist with an unnecessary duplicate in memory.
    del terms
    del generator
    gc.collect()

    if verify:
        if progress:
            print("verifying native cache by loading it back", file=sys.stderr)

        def verify_progress(current: int, total: int) -> None:
            if not progress:
                return
            percentage = 100.0 if total == 0 else 100.0 * current / total
            print(
                f"native verify: {current:,}/{total:,} bytes ({percentage:5.1f}%)",
                file=sys.stderr,
                flush=True,
            )

        restored = load_generator(
            definition,
            target,
            progress=verify_progress if progress else None,
        )
        if not isinstance(restored, EnotsWolleyGenerator):
            raise LegacyEWCacheError("native cache restored the wrong generator type")
        if len(restored.terms) != term_count:
            raise LegacyEWCacheError("native cache term count changed during migration")
        if (restored.terms[-1] if restored.terms else None) != last_term:
            raise LegacyEWCacheError("native cache final term changed during migration")
        if restored.terms[: len(EXPECTED_PREFIX)] != EXPECTED_PREFIX[:term_count]:
            raise LegacyEWCacheError("native cache prefix changed during migration")
        if len(restored.used) != term_count:
            raise LegacyEWCacheError("native cache used-set size is inconsistent")
        if restored.smallest_unused != smallest_unused:
            raise LegacyEWCacheError("native cache smallest-unused state changed")
        if restored.radicals is not None:
            raise LegacyEWCacheError("migrated native cache unexpectedly persisted radicals")

    if progress:
        print(
            f"migration complete: {term_count:,} A336957 terms -> {target}",
            file=sys.stderr,
        )
    return target
