import pickle

import pytest

from lex_earliest_seqs import registry
from lex_earliest_seqs.cache import open_run
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


@pytest.mark.parametrize(
    ("family", "indices"),
    [
        ("square", [1, 4, 9, 16, 25, 36]),
        ("power_of_two", [1, 2, 4, 8, 16, 32]),
        ("self_power", [1, 4, 27, 256]),
    ],
)
def test_prime_index_classifiers(family, indices):
    expected = set(indices)
    for index in range(1, 300):
        assert is_retained_prime_index(index, family) == (index in expected)


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

    for value in (2, 4, 7, 8, 14, 23, 28, 49, 161, 529):
        assert square.allows(value)
    for value in (3, 5, 6, 10, 21, 35, 46):
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
        ("X000014", "prime-coordinate-self-power-ew", "self_power", 12),
    ],
)
def test_registered_sparse_prime_index_sequences(sequence_id, alias, family, count):
    definition = registry.resolve(sequence_id)
    assert definition.oeis is None
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
    ("family", "count"),
    [
        ("square", 20),
        ("power_of_two", 20),
        # p_{i^i} reaches p_27=103 almost immediately; the direct integer
        # scanner is intentionally only a tiny oracle for this family.
        ("self_power", 4),
    ],
)
def test_simple_generator_matches_direct_scanner(family, count):
    generator = SparsePrimeIndexOnlyEnotsWolleyGenerator(family=family)
    reference = ReferenceSparsePrimeIndexOnlyEnotsWolleyGenerator(family=family)

    generator.extend_to(count)
    reference.extend_to(count)

    assert generator.terms == reference.terms
    assert generator.used == set(generator.terms)


@pytest.mark.parametrize(
    ("family", "first_count", "second_count"),
    [
        ("square", 50, 80),
        ("power_of_two", 50, 80),
        ("self_power", 8, 12),
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


def test_registered_definition_constants_have_expected_ids():
    assert SQUARE_INDEX_PRIME_ONLY_ENOTS_WOLLEY.id == "X000012"
    assert POWER_OF_TWO_INDEX_PRIME_ONLY_ENOTS_WOLLEY.id == "X000013"
    assert SELF_POWER_INDEX_PRIME_ONLY_ENOTS_WOLLEY.id == "X000014"


def test_unknown_family_rejected():
    with pytest.raises(ValueError):
        SparsePrimeIndexOnlyEnotsWolleyGenerator(family="not-a-family")
