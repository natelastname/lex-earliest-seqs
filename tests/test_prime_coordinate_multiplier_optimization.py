import pickle

import pytest

from lex_earliest_seqs import registry
from lex_earliest_seqs.zoo.enots_wolley import prime_support
from lex_earliest_seqs.zoo.every_kth_prime_only_enots_wolley import (
    EveryKthPrimeOnlyEnotsWolleyGenerator,
    EveryKthPrimeOnlyPolicy,
)


class _SingleTableBaselineGenerator(EveryKthPrimeOnlyEnotsWolleyGenerator):
    """Version-2 behavior: never switch to the odd multiplier view."""

    def _multiplier_table_for_forbidden(self, forbidden_radical):
        return self.multiplier_values, self.multiplier_successors


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
        29,
        32,
        34,
        41,
        49,
        53,
        56,
        58,
        64,
        67,
        68,
    ]


def test_k_three_odd_multiplier_table_is_exact_odd_submonoid_prefix():
    generator = EveryKthPrimeOnlyEnotsWolleyGenerator(k=3)
    generator._extend_multiplier_table(100)

    expected = [
        value
        for value in generator.multiplier_values
        if value <= 100 and value % 2 == 1
    ]
    actual = [value for value in generator.odd_multiplier_values if value <= 100]

    assert actual == expected
    assert actual == [1, 7, 17, 29, 41, 49, 53, 67, 79, 97]


def test_even_local_forbidden_radical_selects_odd_multiplier_table():
    generator = EveryKthPrimeOnlyEnotsWolleyGenerator(k=4)

    full_values, full_successors = generator._multiplier_table_for_forbidden(11)
    odd_values, odd_successors = generator._multiplier_table_for_forbidden(22)

    assert full_values is generator.multiplier_values
    assert full_successors is generator.multiplier_successors
    assert odd_values is generator.odd_multiplier_values
    assert odd_successors is generator.odd_multiplier_successors


def test_multiplier_tables_never_contain_forbidden_prime_factors():
    for k in (2, 3, 4, 5):
        generator = EveryKthPrimeOnlyEnotsWolleyGenerator(k=k)
        generator._extend_multiplier_table(5_000)
        policy = EveryKthPrimeOnlyPolicy(k)

        assert generator.multiplier_values[0] == 1
        assert generator.odd_multiplier_values[0] == 1
        assert all(policy.allows(value) for value in generator.multiplier_values[1:])
        assert all(value & 1 for value in generator.odd_multiplier_values)
        assert all(
            policy.allows(value) for value in generator.odd_multiplier_values[1:]
        )


def test_multiplier_tables_are_shared_across_stream_primes():
    generator = EveryKthPrimeOnlyEnotsWolleyGenerator(k=4)
    generator.extend_to(200)

    assert len(generator.multiplier_values) > 1
    assert len(generator.odd_multiplier_values) > 1
    assert len(generator.multiplier_successors) > 1
    assert generator.odd_multiplier_successors

    for maps in (generator.multiplier_successors, generator.odd_multiplier_successors):
        for parents in maps.values():
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
def test_odd_fast_path_matches_single_table_streams(k):
    optimized = EveryKthPrimeOnlyEnotsWolleyGenerator(k=k)
    baseline = _SingleTableBaselineGenerator(k=k)

    optimized.extend_to(2_000)
    baseline.extend_to(2_000)

    assert optimized.terms == baseline.terms


def test_registered_prime_coordinate_family_uses_generator_version_three():
    for sequence_id in ("X000009", "X000010", "X000011"):
        assert registry.resolve(sequence_id).generator_version == 3


def test_optimized_generator_pickle_preserves_both_multiplier_views():
    generator = EveryKthPrimeOnlyEnotsWolleyGenerator(k=4)
    generator.extend_to(500)

    values = list(generator.multiplier_values)
    odd_values = list(generator.odd_multiplier_values)
    limit = generator.multiplier_limit
    successors = {
        prime: dict(parents)
        for prime, parents in generator.multiplier_successors.items()
    }
    odd_successors = {
        prime: dict(parents)
        for prime, parents in generator.odd_multiplier_successors.items()
    }

    restored = pickle.loads(pickle.dumps(generator))
    assert restored.k == 4
    assert restored.multiplier_values == values
    assert restored.odd_multiplier_values == odd_values
    assert restored.multiplier_limit == limit
    assert restored.multiplier_successors == successors
    assert restored.odd_multiplier_successors == odd_successors

    restored.extend_to(1_000)
    fresh = EveryKthPrimeOnlyEnotsWolleyGenerator(k=4)
    fresh.extend_to(1_000)
    assert restored.terms == fresh.terms
