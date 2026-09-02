import pickle

from lex_earliest_seqs import registry
from lex_earliest_seqs.cache import open_run
from lex_earliest_seqs.zoo.even_index_prime_enots_wolley import (
    EvenIndexPrimeEnotsWolleyGenerator,
    EvenIndexPrimePolicy,
    ReferenceEvenIndexPrimeEnotsWolleyGenerator,
    is_even_index_prime,
)


def test_even_index_prime_classifier_uses_one_based_prime_indices():
    expected = {
        2: False,
        3: True,
        5: False,
        7: True,
        11: False,
        13: True,
        17: False,
        19: True,
        23: False,
        29: True,
    }
    assert {prime: is_even_index_prime(prime) for prime in expected} == expected
    assert not is_even_index_prime(1)
    assert not is_even_index_prime(9)


def test_even_index_prime_policy_allows_odd_indexed_primes_as_cofactors():
    policy = EvenIndexPrimePolicy()

    for value in (3, 6, 7, 15, 26, 35, 49, 65):
        assert policy.allows(value)

    for value in (1, 2, 5, 10, 11, 22, 25, 55):
        assert not policy.allows(value)


def test_even_index_prime_enots_wolley_reference_prefix():
    definition = registry.resolve("X000006")
    assert definition.oeis is None
    assert registry.resolve("alternating-prime-ew") is definition
    assert "prime-exponents" in definition.projections

    run = open_run(definition, use_cache=False)
    run.ensure(50)
    assert list(run.terms) == [
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
        58,
        145,
        105,
        36,
        74,
        185,
        135,
        42,
        76,
        209,
        99,
        48,
        86,
        215,
        165,
        54,
        98,
        161,
        69,
        60,
    ]


def test_optimized_generator_matches_direct_definition():
    optimized = EvenIndexPrimeEnotsWolleyGenerator()
    reference = ReferenceEvenIndexPrimeEnotsWolleyGenerator()

    optimized.extend_to(2_000)
    reference.extend_to(2_000)

    assert optimized.terms == reference.terms
    assert optimized.used == set(optimized.terms)
    assert optimized.multiplier_successors


def test_even_index_prime_generator_pickle_resumes():
    generator = EvenIndexPrimeEnotsWolleyGenerator()
    generator.extend_to(500)

    restored = pickle.loads(pickle.dumps(generator))
    restored.extend_to(1_000)

    reference = ReferenceEvenIndexPrimeEnotsWolleyGenerator()
    reference.extend_to(1_000)
    assert restored.terms == reference.terms
