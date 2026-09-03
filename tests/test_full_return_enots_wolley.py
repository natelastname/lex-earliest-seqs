import pickle

import pytest

from lex_earliest_seqs import registry
from lex_earliest_seqs.cache import open_run
from lex_earliest_seqs.zoo.full_return_enots_wolley import (
    FULL_RETURN_EW_2_3,
    FULL_RETURN_EW_2_5,
    FULL_RETURN_EW_3_5,
    FullReturnEnotsWolleyGenerator,
    ReferenceFullReturnEnotsWolleyGenerator,
    _multiplier_introduces_new_prime,
    full_return_candidate_allowed,
    full_return_restriction_active,
    make_full_return_enots_wolley_definition,
    target_free,
)


def test_full_return_predicate_uses_only_immediately_previous_term():
    p, q = 2, 3
    previous = 35
    assert target_free(previous, p, q)
    assert full_return_restriction_active(previous, p, q)

    assert full_return_candidate_allowed(55, previous, p, q)
    assert full_return_candidate_allowed(30, previous, p, q)
    assert not full_return_candidate_allowed(14, previous, p, q)
    assert not full_return_candidate_allowed(21, previous, p, q)

    assert not full_return_restriction_active(30, p, q)
    assert full_return_candidate_allowed(14, 30, p, q)
    assert full_return_candidate_allowed(21, 30, p, q)


def test_multiplier_new_prime_test_is_exact_for_stream_cofactors():
    previous_primes = frozenset({2, 3, 7})
    for multiplier in range(1, 500):
        stripped = multiplier
        for prime in previous_primes:
            while stripped % prime == 0:
                stripped //= prime
        assert _multiplier_introduces_new_prime(multiplier, previous_primes) == (
            stripped > 1
        )


@pytest.mark.parametrize(
    ("sequence_id", "alias", "pair", "expected"),
    [
        (
            "X000012",
            "fr-ew-2-3",
            (2, 3),
            [
                1, 2, 6, 15, 35, 77, 66, 10, 65, 91,
                42, 20, 55, 143, 78, 14, 119, 85, 30, 21,
                133, 95, 60, 22, 187, 221, 156, 28, 161, 115,
                90, 26, 247, 209, 132, 34, 323, 437, 138, 33,
                275, 145, 174, 38, 475, 155, 186, 39, 299, 253,
            ],
        ),
        (
            "X000013",
            "fr-ew-2-5",
            (2, 5),
            [
                1, 2, 6, 15, 35, 14, 12, 33, 77, 70,
                18, 39, 91, 119, 51, 30, 22, 143, 117, 21,
                133, 190, 24, 63, 161, 230, 26, 221, 153, 57,
                209, 110, 28, 147, 69, 253, 187, 170, 36, 87,
                203, 140, 34, 323, 171, 60, 44, 319, 261, 90,
            ],
        ),
        (
            "X000014",
            "fr-ew-3-5",
            (3, 5),
            [
                1, 2, 14, 77, 143, 26, 28, 105, 33, 22,
                34, 119, 91, 52, 30, 21, 133, 38, 44, 165,
                35, 56, 46, 253, 187, 68, 58, 203, 161, 92,
                60, 39, 221, 136, 62, 217, 259, 74, 76, 209,
                319, 116, 82, 287, 301, 86, 88, 341, 403, 104,
            ],
        ),
    ],
)
def test_registered_full_return_prefixes(sequence_id, alias, pair, expected):
    definition = registry.resolve(sequence_id)
    assert definition.oeis is None
    assert definition.definition_version == 2
    assert registry.resolve(alias) is definition
    assert "prime-exponents" in definition.projections

    generator = definition.generator_factory()
    assert isinstance(generator, FullReturnEnotsWolleyGenerator)
    assert (generator.p, generator.q) == pair

    run = open_run(definition, use_cache=False)
    run.ensure(len(expected))
    assert list(run.terms) == expected


def test_registered_definition_constants_have_expected_ids():
    assert FULL_RETURN_EW_2_3.id == "X000012"
    assert FULL_RETURN_EW_2_5.id == "X000013"
    assert FULL_RETURN_EW_3_5.id == "X000014"


@pytest.mark.parametrize("pair", [(2, 3), (2, 5), (3, 5), (2, 7), (5, 7)])
def test_optimized_generator_matches_direct_definition(pair):
    p, q = pair
    optimized = FullReturnEnotsWolleyGenerator(p=p, q=q)
    reference = ReferenceFullReturnEnotsWolleyGenerator(p=p, q=q)

    optimized.extend_to(5_000)
    reference.extend_to(5_000)

    assert optimized.terms == reference.terms
    assert optimized.used == set(optimized.terms)
    assert optimized.unused_multiplier_successors
    # The specialized hot path derives radicals directly from prime supports;
    # it must not allocate ordinary EW's dense candidate-value radical table.
    assert optimized.radicals is None


@pytest.mark.parametrize("pair", [(2, 3), (2, 5), (3, 5)])
def test_every_target_return_after_one_free_term_is_full(pair):
    p, q = pair
    generator = FullReturnEnotsWolleyGenerator(p=p, q=q)
    generator.extend_to(5_000)

    for index in range(2, len(generator.terms)):
        previous = generator.terms[index - 1]
        value = generator.terms[index]
        if full_return_restriction_active(previous, p, q):
            has_p = value % p == 0
            has_q = value % q == 0
            assert has_p == has_q


def test_state_local_rejections_can_become_eligible_later():
    cases = [
        ((2, 3), 14, 16),
        ((2, 5), 55, 71),
        ((3, 5), 6, 58),
    ]
    for (p, q), value, expected_position in cases:
        generator = FullReturnEnotsWolleyGenerator(p=p, q=q)
        generator.extend_to(expected_position)
        assert generator.terms[expected_position - 1] == value


def test_full_return_generator_pickle_resumes():
    generator = FullReturnEnotsWolleyGenerator(p=3, q=5)
    generator.extend_to(500)
    restored = pickle.loads(pickle.dumps(generator))

    assert (restored.p, restored.q) == (3, 5)
    restored.extend_to(1_000)

    reference = ReferenceFullReturnEnotsWolleyGenerator(p=3, q=5)
    reference.extend_to(1_000)
    assert restored.terms == reference.terms


def test_pair_order_is_canonicalized():
    generator = FullReturnEnotsWolleyGenerator(p=5, q=2)
    assert (generator.p, generator.q) == (2, 5)

    definition = make_full_return_enots_wolley_definition(
        id="X999999",
        p=5,
        q=2,
    )
    instance = definition.generator_factory()
    assert (instance.p, instance.q) == (2, 5)


@pytest.mark.parametrize("pair", [(2, 2), (4, 5), (3, 9)])
def test_invalid_target_pairs_are_rejected(pair):
    with pytest.raises(ValueError):
        FullReturnEnotsWolleyGenerator(p=pair[0], q=pair[1])


def test_noninteger_target_pair_is_rejected():
    with pytest.raises(TypeError):
        FullReturnEnotsWolleyGenerator(p=2.0, q=3)
