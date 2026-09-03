import pickle

import pytest

from lex_earliest_seqs import registry
from lex_earliest_seqs.zoo.enots_wolley import prime_support
from lex_earliest_seqs.zoo.every_kth_prime_only_enots_wolley import (
    EveryKthPrimeOnlyEnotsWolleyGenerator,
    EveryKthPrimeOnlyPolicy,
    ReferenceEveryKthPrimeOnlyEnotsWolleyGenerator,
)


def test_k_three_multiplier_table_is_exact_retained_prime_monoid_prefix():
    generator = EveryKthPrimeOnlyEnotsWolleyGenerator(k=3)
    generator._extend_multiplier_table(100)

    expected = [
        value
        for value in range(1, 101)
        if value == 1 or EveryKthPrimeOnlyPolicy(3).allows(value)
    ]
    actual = [value for value in generator.multiplier_values if value <= 100]
    assert actual == expected
    assert actual[:20] == [
        1,
        2,
        4,
        7,
        8,
        14,
        16,
        17,
        28,
        32,
        34,
        49,
        56,
        64,
        68,
        98,
    ]


def test_multiplier_table_never_contains_forbidden_prime_factors():
    for k in (2, 3, 4, 5):
        generator = EveryKthPrimeOnlyEnotsWolleyGenerator(k=k)
        generator._extend_multiplier_table(5_000)
        policy = EveryKthPrimeOnlyPolicy(k)
        assert generator.multiplier_values[0] == 1
        assert all(policy.allows(value) for value in generator.multiplier_values[1:])


def test_multiplier_table_is_shared_across_stream_primes():
    generator = EveryKthPrimeOnlyEnotsWolleyGenerator(k=4)
    generator.extend_to(200)

    assert len(generator.multiplier_values) > 1
    assert len(generator.multiplier_successors) > 1
    # Successor maps store indices into one common multiplier_values list.
    for parents in generator.multiplier_successors.values():
        assert all(type(index) is int for index in parents)
        assert all(type(successor) is int for successor in parents.values())
        assert all(successor > index for index, successor in parents.items())


@pytest.mark.parametrize("k", [1, 2, 3, 4, 5])
def test_stripped_new_prime_test_matches_support_definition(k):
    generator = EveryKthPrimeOnlyEnotsWolleyGenerator(k=k)
    generator._extend_multiplier_table(2_000)

    predecessor_supports = [
        frozenset({2}),
        frozenset({2, 7}),
        frozenset({2, 7, 17}),
    ]
    for previous_support in predecessor_supports:
        for multiplier in generator.multiplier_values[:500]:
            expected = bool(prime_support(multiplier) - previous_support)
            assert (
                generator._multiplier_introduces_new_prime(
                    multiplier,
                    previous_support,
                )
                == expected
            )


@pytest.mark.parametrize("k", [2, 3, 4, 5])
def test_retained_multiplier_generator_matches_direct_definition_at_2000_terms(k):
    optimized = EveryKthPrimeOnlyEnotsWolleyGenerator(k=k)
    reference = ReferenceEveryKthPrimeOnlyEnotsWolleyGenerator(k=k)

    optimized.extend_to(2_000)
    reference.extend_to(2_000)

    assert optimized.terms == reference.terms
    assert optimized.used == set(optimized.terms)


def test_registered_prime_coordinate_family_uses_generator_version_two():
    for sequence_id in ("X000009", "X000010", "X000011"):
        assert registry.resolve(sequence_id).generator_version == 2


def test_optimized_generator_pickle_preserves_multiplier_coordinate_state():
    generator = EveryKthPrimeOnlyEnotsWolleyGenerator(k=4)
    generator.extend_to(500)

    values = list(generator.multiplier_values)
    limit = generator.multiplier_limit
    successors = {
        prime: dict(parents)
        for prime, parents in generator.multiplier_successors.items()
    }

    restored = pickle.loads(pickle.dumps(generator))
    assert restored.k == 4
    assert restored.multiplier_values == values
    assert restored.multiplier_limit == limit
    assert restored.multiplier_successors == successors

    restored.extend_to(1_000)
    reference = ReferenceEveryKthPrimeOnlyEnotsWolleyGenerator(k=4)
    reference.extend_to(1_000)
    assert restored.terms == reference.terms
