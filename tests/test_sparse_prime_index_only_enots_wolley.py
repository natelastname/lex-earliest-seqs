import pickle

import pytest

from lex_earliest_seqs import registry
from lex_earliest_seqs.cache import open_run
from lex_earliest_seqs.zoo import sparse_prime_index_only_enots_wolley as sparse_module
from lex_earliest_seqs.zoo.enots_wolley import prime_support
from lex_earliest_seqs.zoo.sparse_prime_index_only_enots_wolley import (
    POWER_OF_TWO_INDEX_PRIME_ONLY_ENOTS_WOLLEY,
    SELF_POWER_INDEX_PRIME_ONLY_ENOTS_WOLLEY,
    SQUARE_INDEX_PRIME_ONLY_ENOTS_WOLLEY,
    ReferenceSparsePrimeIndexOnlyEnotsWolleyGenerator,
    SparsePrimeIndexOnlyEnotsWolleyGenerator,
    SparsePrimeIndexOnlyPolicy,
    is_retained_prime,
    is_retained_prime_index,
)


def test_prime_index_classifiers():
    for index in range(1, 300):
        assert is_retained_prime_index(index, "square") == (
            int(index**0.5) ** 2 == index
        )
        assert is_retained_prime_index(index, "power_of_two") == (
            index & (index - 1) == 0
        )
        assert is_retained_prime_index(index, "self_power") == (
            index in {1, 4, 27, 256}
        )


def test_retained_prime_prefixes():
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53]

    assert [p for p in primes if is_retained_prime(p, "square")] == [2, 7, 23, 53]
    assert [p for p in primes if is_retained_prime(p, "power_of_two")] == [
        2,
        3,
        7,
        19,
        53,
    ]
    assert [p for p in primes if is_retained_prime(p, "self_power")] == [2, 7]


def test_policies_forbid_nonretained_prime_cofactors():
    square = SparsePrimeIndexOnlyPolicy("square")
    power_two = SparsePrimeIndexOnlyPolicy("power_of_two")
    self_power = SparsePrimeIndexOnlyPolicy("self_power")

    for value in (2, 4, 7, 8, 14, 23, 28, 46, 49, 161, 529):
        assert square.allows(value)
    for value in (3, 5, 6, 10, 21, 22, 35):
        assert not square.allows(value)

    for value in (2, 3, 4, 6, 7, 9, 14, 19, 21, 53):
        assert power_two.allows(value)
    for value in (5, 10, 11, 22, 23, 46):
        assert not power_two.allows(value)

    for value in (2, 4, 7, 8, 14, 49, 103, 206, 721):
        assert self_power.allows(value)
    for value in (3, 5, 6, 11, 22, 23, 46):
        assert not self_power.allows(value)


@pytest.mark.parametrize(
    ("sequence_id", "alias", "family", "count"),
    [
        ("X000012", "prime-coordinate-square-ew", "square", 60),
        ("X000013", "prime-coordinate-power-of-two-ew", "power_of_two", 60),
        ("X000014", "prime-coordinate-self-power-ew", "self_power", 20),
    ],
)
def test_registered_sparse_prime_index_sequences(sequence_id, alias, family, count):
    definition = registry.resolve(sequence_id)
    assert definition.oeis is None
    assert definition.generator_version == 2
    assert registry.resolve(alias) is definition
    assert "prime-exponents" in definition.projections

    generator = definition.generator_factory()
    assert isinstance(generator, SparsePrimeIndexOnlyEnotsWolleyGenerator)
    assert generator.family == family
    assert generator.terms == [1, 2]

    run = open_run(definition, use_cache=False)
    run.ensure(count)
    assert len(run.terms) == count
    assert len(set(run.terms)) == count

    policy = SparsePrimeIndexOnlyPolicy(family)
    for term in run.terms[1:]:
        assert policy.allows(term)


@pytest.mark.parametrize(
    ("family", "expected"),
    [
        (
            "square",
            [
                1,
                2,
                14,
                161,
                1219,
                106,
                28,
                679,
                2231,
                46,
                56,
                371,
                5141,
                194,
                92,
                1127,
                1057,
                302,
                184,
                3703,
            ],
        ),
        (
            "power_of_two",
            [
                1,
                2,
                6,
                21,
                133,
                38,
                12,
                63,
                371,
                106,
                18,
                57,
                931,
                14,
                24,
                159,
                1007,
                76,
                28,
                147,
            ],
        ),
        (
            "self_power",
            [
                1,
                2,
                14,
                721,
                166757,
                3238,
                28,
                5047,
                2954761,
                57374,
                56,
                11333,
                17175971,
                206,
                98,
                79331,
                46444253,
                114748,
                112,
                35329,
            ],
        ),
    ],
)
def test_optimized_generator_preserves_known_prefixes(family, expected):
    generator = SparsePrimeIndexOnlyEnotsWolleyGenerator(family=family)
    generator.extend_to(len(expected))
    assert generator.terms == expected


@pytest.mark.parametrize(
    ("family", "count"),
    [
        ("square", 20),
        ("power_of_two", 20),
        ("self_power", 4),
    ],
)
def test_optimized_generator_matches_direct_scanner(family, count):
    generator = SparsePrimeIndexOnlyEnotsWolleyGenerator(family=family)
    reference = ReferenceSparsePrimeIndexOnlyEnotsWolleyGenerator(family=family)

    generator.extend_to(count)
    reference.extend_to(count)

    assert generator.terms == reference.terms
    assert generator.used == set(generator.terms)


@pytest.mark.parametrize(
    ("family", "count"),
    [
        ("square", 80),
        ("power_of_two", 80),
        ("self_power", 20),
    ],
)
def test_support_masks_match_factorization(family, count):
    generator = SparsePrimeIndexOnlyEnotsWolleyGenerator(family=family)
    generator.extend_to(count)

    assert len(generator.term_support_masks) == len(generator.terms)
    for term, mask in zip(generator.terms, generator.term_support_masks, strict=True):
        expected = {
            generator.retained_primes[position]
            for position in range(mask.bit_length())
            if mask & (1 << position)
        }
        assert expected == set(prime_support(term))

    for value, mask in zip(
        generator.multiplier_values,
        generator.multiplier_support_masks,
        strict=True,
    ):
        expected = {
            generator.retained_primes[position]
            for position in range(mask.bit_length())
            if mask & (1 << position)
        }
        assert expected == set(prime_support(value))


@pytest.mark.parametrize(
    ("family", "count"),
    [
        ("square", 100),
        ("power_of_two", 100),
        ("self_power", 20),
    ],
)
def test_production_candidate_engine_does_not_factor_terms(monkeypatch, family, count):
    def fail_prime_support(_value):
        raise AssertionError("production sparse generator unexpectedly factored a term")

    generator = SparsePrimeIndexOnlyEnotsWolleyGenerator(family=family)
    monkeypatch.setattr(sparse_module, "prime_support", fail_prime_support)
    generator.extend_to(count)
    assert len(generator.terms) == count


@pytest.mark.parametrize(
    ("family", "first_count", "second_count"),
    [
        ("square", 50, 80),
        ("power_of_two", 50, 80),
        ("self_power", 12, 20),
    ],
)
def test_generator_pickle_resumes(family, first_count, second_count):
    generator = SparsePrimeIndexOnlyEnotsWolleyGenerator(family=family)
    generator.extend_to(first_count)

    restored = pickle.loads(pickle.dumps(generator))
    assert restored.family == family
    assert restored.policy == SparsePrimeIndexOnlyPolicy(family)
    restored.extend_to(second_count)

    fresh = SparsePrimeIndexOnlyEnotsWolleyGenerator(family=family)
    fresh.extend_to(second_count)
    assert restored.terms == fresh.terms
    assert restored.term_support_masks == fresh.term_support_masks


def test_pickle_resume_reseeds_high_index_prime_cursor(monkeypatch):
    generator = SparsePrimeIndexOnlyEnotsWolleyGenerator(family="self_power")
    generator.extend_to(20)
    restored = pickle.loads(pickle.dumps(generator))

    checkpoint_position = len(restored.retained_primes) - 1
    checkpoint = (
        restored._allowed_prime_index(checkpoint_position),
        restored.retained_primes[checkpoint_position],
    )
    next_position = len(restored.retained_primes)
    calls: list[tuple[int, int]] = []

    monkeypatch.setattr(
        sparse_module,
        "remember_nth_prime",
        lambda index, prime: calls.append((index, prime)),
    )
    monkeypatch.setattr(sparse_module, "nth_prime", lambda _index: 999_983)

    assert restored._prime_at_position(next_position) == 999_983
    assert calls == [checkpoint]


def test_registered_definition_constants_have_expected_ids():
    assert SQUARE_INDEX_PRIME_ONLY_ENOTS_WOLLEY.id == "X000012"
    assert POWER_OF_TWO_INDEX_PRIME_ONLY_ENOTS_WOLLEY.id == "X000013"
    assert SELF_POWER_INDEX_PRIME_ONLY_ENOTS_WOLLEY.id == "X000014"


def test_unknown_family_rejected():
    with pytest.raises(ValueError):
        SparsePrimeIndexOnlyEnotsWolleyGenerator(family="not-a-family")
