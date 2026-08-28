from __future__ import annotations

import pickle

import pytest

from lex_earliest_seqs import registry
from lex_earliest_seqs._legacy_ew_migration import (
    LEGACY_FORMAT_VERSION,
    LEGACY_SEQUENCE_ID,
    LegacyEWCacheError,
    migrate_legacy_ew_cache,
)
from lex_earliest_seqs.cache import load_generator
from lex_earliest_seqs.zoo.enots_wolley import EnotsWolleyGenerator


def _write_legacy_cache(path, terms, *, sequence_id=LEGACY_SEQUENCE_ID):
    with path.open("wb") as handle:
        pickle.dump(
            {
                "format_version": LEGACY_FORMAT_VERSION,
                "sequence_id": sequence_id,
                "term_count": len(terms),
                "terms": terms,
            },
            handle,
            protocol=pickle.HIGHEST_PROTOCOL,
        )


def test_legacy_cache_migrates_without_replaying_and_continues(tmp_path):
    expected = EnotsWolleyGenerator()
    expected.extend_to(600)
    legacy_terms = expected.terms[:300]

    source = tmp_path / "terms-v1.pkl"
    destination = tmp_path / "A336957.pkl"
    _write_legacy_cache(source, legacy_terms)

    migrate_legacy_ew_cache(
        source,
        destination,
        progress=False,
    )

    definition = registry.resolve("A336957")
    restored = load_generator(definition, destination)
    assert isinstance(restored, EnotsWolleyGenerator)
    assert restored.terms == legacy_terms
    assert restored.used == set(legacy_terms)
    assert restored.radicals is None
    assert restored.smallest_unused not in restored.used
    assert all(value in restored.used for value in range(1, restored.smallest_unused))

    restored.extend_to(600)
    assert restored.terms == expected.terms


def test_legacy_cache_rejects_wrong_sequence_identity(tmp_path):
    source = tmp_path / "wrong.pkl"
    _write_legacy_cache(source, [1, 2], sequence_id="not-a336957")

    with pytest.raises(LegacyEWCacheError, match="sequence_id"):
        migrate_legacy_ew_cache(
            source,
            tmp_path / "A336957.pkl",
            progress=False,
        )
