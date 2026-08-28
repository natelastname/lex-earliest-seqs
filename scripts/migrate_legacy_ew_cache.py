#!/usr/bin/env python3
"""One-off migration of the historical enots-wolley-2 A336957 cache.

Run from the repository root with:

    uv run python scripts/migrate_legacy_ew_cache.py

Delete this script and ``lex_earliest_seqs._legacy_ew_migration`` after the
large legacy cache has been converted and verified.
"""

from __future__ import annotations

from pathlib import Path

from cyclopts import run

from lex_earliest_seqs._legacy_ew_migration import (
    DEFAULT_LEGACY_CACHE_PATH,
    migrate_legacy_ew_cache,
)


def main(
    source: Path = DEFAULT_LEGACY_CACHE_PATH,
    destination: Path | None = None,
    *,
    force: bool = False,
    verify: bool = True,
    progress: bool = True,
) -> None:
    """Convert the old A336957 term-list pickle into the native generator pickle.

    Parameters
    ----------
    source
        Historical enots-wolley-2 cache. Defaults to
        ``~/.cache/enots-wolley-2/terms-v1.pkl``.
    destination
        Native cache path. Defaults to the normal A336957 cache path under
        ``~/.cache/lex-earliest-seqs`` (or ``$XDG_CACHE_HOME``).
    force
        Replace an existing native A336957 cache.
    verify
        Load the newly written native pickle back and verify reconstructed state.
    progress
        Print migration progress to stderr.
    """

    migrate_legacy_ew_cache(
        source=source,
        destination=destination,
        force=force,
        verify=verify,
        progress=progress,
    )


if __name__ == "__main__":
    run(main)
