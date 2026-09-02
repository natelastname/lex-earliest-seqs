import pickle

import pytest

from lex_earliest_seqs import registry
from lex_earliest_seqs.cache import open_run
from lex_earliest_seqs.zoo.enots_wolley import EnotsWolleyGenerator
from lex_earliest_seqs.zoo.every_kth_prime_enots_wolley import (
    EVERY_FOURTH_PRIME_ENOTS_WOLLEY,
    EVERY_SECOND_PRIME_ENOTS_WOLLEY,
    EVERY_THIRD_PRIME_ENOTS_WOLLEY,
    EveryKthPrimeEnotsWolleyGenerator,
    EveryKthPrimePolicy,
    ReferenceEveryKthPrimeEnotsWolleyGenerator,
    is_every_kth_prime,
    nth_prime,
    prime_index,
)


def test_prime_index_helpers_use_one_based_indices():
    expected = {
        2: 1,
        3: 2,
        5: 3,
        7: 4,
        11: 5,
        13: 6,
        17: 7,
        19: 8,
        23: 9,
        29: 10,
    }
    assert {prime: prime_index(prime) for prime in expected} == expected
    assert prime_index(1) is None
    assert prime_index(9) is None
    assert [nth_prime(index) for index in range(1, 11)] == list(expected)


def test_every_kth_prime_classifier_is_generic():
    assert [prime for prime in (2, 3, 5, 7, 11, 13, 17, 19) if is_every_kth_prime(prime, 2)] == [
        3,
        7,
        13,
        19,
    ]
    assert [prime for prime in (2, 3, 5, 7, 11, 13, 17, 19, 23) if is_every_kth_prime(prime, 3)] == [
        5,
        13,
        23,
    ]
    assert is_every_kth_prime(11, 5)
    assert not is_every_kth_prime(7, 5)
    assert not is_every_kth_prime(25, 5)


def test_every_kth_prime_policy_keeps_other_primes_as_cofactors():
    policy = EveryKthPrimePolicy(3)

    for value in (5, 10, 13, 15, 26, 39, 65):
        assert policy.allows(value)

    for value in (1, 2, 3, 6, 7, 11, 21, 22, 49, 77):
        assert not policy.allows(value)


@pytest.mark.parametrize(
    ("sequence_id", "alias", "k", "expected"),
    [
        (
            "X000006",
            "alternating-prime-ew",
            2,
            [
                1,
                3,
                6,
                14,
                35,
                15,
                12,
                26,
                65,
                45,
                18,
                28,
                77,
                33,
                24,
                38,
                95,
                75,
                21,
                56,
                52,
                39,
                51,
                119,
                70,
                30,
                57,
                133,
                91,
                78,
            ],
        ),
        (
            "X000007",
            "every-third-prime-ew",
            3,
            [
                1,
                5,
                10,
                26,
                39,
                15,
                20,
                46,
                69,
                45,
                35,
                91,
                52,
                30,
                55,
                143,
                78,
                40,
                85,
                221,
                104,
                50,
                75,
                111,
                74,
                70,
                65,
                117,
                138,
                80,
            ],
        ),
        (
            "X000008",
            "every-fourth-prime-ew",
            4,
            [
                1,
                7,
                14,
                38,
                57,
                21,
                28,
                74,
                111,
                63,
                35,
                95,
                76,
                42,
                77,
                209,
                114,
                56,
                91,
                247,
                152,
                70,
                105,
                159,
                106,
                98,
                119,
                323,
                171,
                84,
            ],
        ),
    ],
)
def test_registered_every_kth_prime_prefixes(sequence_id, alias, k, expected):
    definition = registry.resolve(sequence_id)
    assert definition.oeis is None
    assert registry.resolve(alias) is definition
    assert "prime-exponents" in definition.projections

    generator = definition.generator_factory()
    assert isinstance(generator, EveryKthPrimeEnotsWolleyGenerator)
    assert generator.k == k
    assert generator.terms == [1, nth_prime(k)]

    run = open_run(definition, use_cache=False)
    run.ensure(len(expected))
    assert list(run.terms) == expected


def test_registered_definition_constants_have_expected_ids():
    assert EVERY_SECOND_PRIME_ENOTS_WOLLEY.id == "X000006"
    assert EVERY_THIRD_PRIME_ENOTS_WOLLEY.id == "X000007"
    assert EVERY_FOURTH_PRIME_ENOTS_WOLLEY.id == "X000008"


@pytest.mark.parametrize("k", [2, 3, 4, 5])
def test_optimized_generator_matches_direct_definition_for_arbitrary_k(k):
    optimized = EveryKthPrimeEnotsWolleyGenerator(k=k)
    reference = ReferenceEveryKthPrimeEnotsWolleyGenerator(k=k)

    optimized.extend_to(500)
    reference.extend_to(500)

    assert optimized.terms == reference.terms
    assert optimized.used == set(optimized.terms)
    assert optimized.multiplier_successors


def test_k_one_reproduces_ordinary_enots_wolley():
    generalized = EveryKthPrimeEnotsWolleyGenerator(k=1)
    ordinary = EnotsWolleyGenerator()

    generalized.extend_to(500)
    ordinary.extend_to(500)

    assert generalized.terms == ordinary.terms


def test_unregistered_k_five_uses_p5_seed_without_special_case():
    generator = EveryKthPrimeEnotsWolleyGenerator(k=5)
    generator.extend_to(20)
    assert generator.terms == [
        1,
        11,
        22,
        58,
        87,
        33,
        44,
        94,
        141,
        99,
        55,
        145,
        116,
        66,
        77,
        203,
        174,
        88,
        143,
        377,
    ]


def test_every_kth_prime_generator_pickle_resumes_with_k():
    generator = EveryKthPrimeEnotsWolleyGenerator(k=4)
    generator.extend_to(250)

    restored = pickle.loads(pickle.dumps(generator))
    assert restored.k == 4
    assert restored.policy == EveryKthPrimePolicy(4)
    restored.extend_to(500)

    reference = ReferenceEveryKthPrimeEnotsWolleyGenerator(k=4)
    reference.extend_to(500)
    assert restored.terms == reference.terms


@pytest.mark.parametrize("bad_k", [0, -1])
def test_every_kth_prime_generator_rejects_nonpositive_k(bad_k):
    with pytest.raises(ValueError):
        EveryKthPrimeEnotsWolleyGenerator(k=bad_k)


def test_every_kth_prime_generator_rejects_noninteger_k():
    with pytest.raises(TypeError):
        EveryKthPrimeEnotsWolleyGenerator(k=2.0)
