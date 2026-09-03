import pickle
from math import isqrt

import pytest

from lex_earliest_seqs import registry
from lex_earliest_seqs.cache import open_run
from lex_earliest_seqs.zoo.enots_wolley import EnotsWolleyGenerator, prime_support
from lex_earliest_seqs.zoo.every_kth_prime_only_enots_wolley import (
    EVERY_FOURTH_PRIME_ONLY_ENOTS_WOLLEY,
    EVERY_SECOND_PRIME_ONLY_ENOTS_WOLLEY,
    EVERY_THIRD_PRIME_ONLY_ENOTS_WOLLEY,
    EveryKthPrimeOnlyEnotsWolleyGenerator,
    EveryKthPrimeOnlyPolicy,
    ReferenceEveryKthPrimeOnlyEnotsWolleyGenerator,
    is_retained_prime,
)


def _is_prime(value: int) -> bool:
    if value < 2:
        return False
    for divisor in range(2, isqrt(value) + 1):
        if value % divisor == 0:
            return False
    return True


def _reference_prime_index(prime: int) -> int:
    assert _is_prime(prime)
    return sum(1 for value in range(2, prime + 1) if _is_prime(value))


def _reference_retained_prime(prime: int, k: int) -> bool:
    return _is_prime(prime) and (_reference_prime_index(prime) - 1) % k == 0


def _reference_allowed(value: int, k: int) -> bool:
    support = prime_support(value)
    return bool(support) and all(_reference_retained_prime(prime, k) for prime in support)


def test_retained_prime_classifier_starts_at_two_and_strides_by_k():
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]

    assert [prime for prime in primes if is_retained_prime(prime, 2)] == [
        2,
        5,
        11,
        17,
        23,
        31,
    ]
    assert [prime for prime in primes if is_retained_prime(prime, 3)] == [
        2,
        7,
        17,
        29,
    ]
    assert [prime for prime in primes if is_retained_prime(prime, 4)] == [
        2,
        11,
        23,
    ]

    for k in (2, 3, 4, 5):
        for value in range(1, 100):
            expected = _reference_retained_prime(value, k) if _is_prime(value) else False
            assert is_retained_prime(value, k) == expected


def test_k_three_policy_forbids_three_and_five_even_as_cofactors():
    policy = EveryKthPrimeOnlyPolicy(3)

    for value in (2, 4, 7, 8, 14, 17, 28, 34, 49, 119, 289):
        assert policy.allows(value)

    for value in (1, 3, 5, 6, 10, 15, 21, 35, 42, 57, 63, 111):
        assert not policy.allows(value)


def test_policy_matches_independent_prime_index_definition():
    for k in (1, 2, 3, 4, 5):
        policy = EveryKthPrimeOnlyPolicy(k)
        for value in range(2, 500):
            assert policy.allows(value) == _reference_allowed(value, k)


@pytest.mark.parametrize(
    ("sequence_id", "alias", "k", "expected"),
    [
        (
            "X000009",
            "prime-coordinate-ew-2",
            2,
            [
                1,
                2,
                10,
                55,
                187,
                34,
                20,
                115,
                253,
                22,
                40,
                85,
                391,
                46,
                44,
                275,
                155,
                62,
                68,
                425,
                205,
                82,
                88,
                341,
                527,
                136,
                50,
                235,
                517,
                176,
            ],
        ),
        (
            "X000010",
            "prime-coordinate-ew-3",
            3,
            [
                1,
                2,
                14,
                119,
                493,
                58,
                28,
                287,
                697,
                34,
                56,
                203,
                1189,
                82,
                68,
                833,
                371,
                106,
                116,
                1421,
                469,
                134,
                136,
                901,
                1537,
                232,
                98,
                553,
                1343,
                272,
            ],
        ),
        (
            "X000011",
            "prime-coordinate-ew-4",
            4,
            [
                1,
                2,
                22,
                253,
                943,
                82,
                44,
                649,
                1357,
                46,
                88,
                451,
                2419,
                118,
                92,
                1679,
                803,
                176,
                164,
                2993,
                4307,
                236,
                184,
                2231,
                1067,
                242,
                146,
                7081,
                3977,
                328,
            ],
        ),
    ],
)
def test_registered_prime_coordinate_prefixes(sequence_id, alias, k, expected):
    definition = registry.resolve(sequence_id)
    assert definition.oeis is None
    assert registry.resolve(alias) is definition
    assert "prime-exponents" in definition.projections

    generator = definition.generator_factory()
    assert isinstance(generator, EveryKthPrimeOnlyEnotsWolleyGenerator)
    assert generator.k == k
    assert generator.terms == [1, 2]

    run = open_run(definition, use_cache=False)
    run.ensure(len(expected))
    assert list(run.terms) == expected

    for term in run.terms[1:]:
        assert _reference_allowed(term, k)


def test_registered_definition_constants_have_expected_ids():
    assert EVERY_SECOND_PRIME_ONLY_ENOTS_WOLLEY.id == "X000009"
    assert EVERY_THIRD_PRIME_ONLY_ENOTS_WOLLEY.id == "X000010"
    assert EVERY_FOURTH_PRIME_ONLY_ENOTS_WOLLEY.id == "X000011"


@pytest.mark.parametrize("k", [2, 3, 4, 5])
def test_optimized_generator_matches_direct_definition_for_arbitrary_k(k):
    optimized = EveryKthPrimeOnlyEnotsWolleyGenerator(k=k)
    reference = ReferenceEveryKthPrimeOnlyEnotsWolleyGenerator(k=k)

    optimized.extend_to(150)
    reference.extend_to(150)

    assert optimized.terms == reference.terms
    assert optimized.used == set(optimized.terms)
    assert optimized.multiplier_successors


def test_k_one_reproduces_ordinary_enots_wolley():
    generalized = EveryKthPrimeOnlyEnotsWolleyGenerator(k=1)
    ordinary = EnotsWolleyGenerator()

    generalized.extend_to(500)
    ordinary.extend_to(500)

    assert generalized.terms == ordinary.terms


def test_unregistered_k_five_is_generic():
    generator = EveryKthPrimeOnlyEnotsWolleyGenerator(k=5)
    reference = ReferenceEveryKthPrimeOnlyEnotsWolleyGenerator(k=5)

    generator.extend_to(100)
    reference.extend_to(100)
    assert generator.terms == reference.terms


def test_generator_pickle_resumes_with_k():
    generator = EveryKthPrimeOnlyEnotsWolleyGenerator(k=3)
    generator.extend_to(150)

    restored = pickle.loads(pickle.dumps(generator))
    assert restored.k == 3
    assert restored.policy == EveryKthPrimeOnlyPolicy(3)
    restored.extend_to(250)

    reference = ReferenceEveryKthPrimeOnlyEnotsWolleyGenerator(k=3)
    reference.extend_to(250)
    assert restored.terms == reference.terms


@pytest.mark.parametrize("bad_k", [0, -1])
def test_generator_rejects_nonpositive_k(bad_k):
    with pytest.raises(ValueError):
        EveryKthPrimeOnlyEnotsWolleyGenerator(k=bad_k)


def test_generator_rejects_noninteger_k():
    with pytest.raises(TypeError):
        EveryKthPrimeOnlyEnotsWolleyGenerator(k=3.0)
