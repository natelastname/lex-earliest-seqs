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
    full_return_candidate_allowed,
    full_return_restriction_active,
    make_full_return_enots_wolley_definition,
    target_free,
)


def test_full_return_predicate_only_suppresses_one_sided_returns_after_two_free_terms():
    p, q = 2, 3
    previous, two_back = 35, 77
    assert target_free(previous, p, q)
    assert target_free(two_back, p, q)
    assert full_return_restriction_active(previous, two_back, p, q)

    assert full_return_candidate_allowed(55, previous, two_back, p, q)
    assert full_return_candidate_allowed(30, previous, two_back, p, q)
    assert not full_return_candidate_allowed(14, previous, two_back, p, q)
    assert not full_return_candidate_allowed(21, previous, two_back, p, q)

    # As soon as one predecessor term contains a target prime, the extra
    # restriction disappears and proper one-sided target terms are legal again.
    assert not full_return_restriction_active(30, previous, p, q)
    assert full_return_candidate_allowed(14, 30, previous, p, q)
    assert full_return_candidate_allowed(21, 30, previous, p, q)


@pytest.mark.parametrize(
    ("sequence_id", "alias", "pair", "expected"),
    [
        (
            "X000012",
            "fr-ew-2-3",
            (2, 3),
            [
                1, 2, 6, 15, 35, 14, 12, 33, 55, 10,
                18, 21, 77, 22, 20, 45, 39, 26, 28, 63,
                51, 34, 38, 57, 69, 46, 40, 65, 91, 42,
                30, 85, 119, 84, 44, 143, 117, 24, 50, 95,
                133, 126, 52, 221, 153, 36, 56, 161, 115, 60,
            ],
        ),
        (
            "X000013",
            "fr-ew-2-5",
            (2, 5),
            [
                1, 2, 6, 15, 35, 14, 12, 33, 55, 10,
                18, 21, 77, 110, 24, 39, 65, 20, 22, 99,
                45, 40, 26, 91, 63, 30, 34, 119, 105, 36,
                38, 95, 75, 42, 44, 143, 117, 51, 170, 28,
                133, 57, 60, 46, 161, 147, 87, 290, 52, 221,
            ],
        ),
        (
            "X000014",
            "fr-ew-3-5",
            (3, 5),
            [
                1, 2, 14, 77, 143, 26, 28, 105, 33, 22,
                10, 15, 21, 56, 20, 45, 39, 52, 34, 119,
                91, 104, 30, 35, 133, 38, 44, 165, 51, 68,
                40, 55, 99, 6, 46, 115, 65, 78, 42, 161,
                253, 88, 58, 203, 217, 62, 60, 57, 209, 110,
            ],
        ),
    ],
)
def test_registered_full_return_prefixes(sequence_id, alias, pair, expected):
    definition = registry.resolve(sequence_id)
    assert definition.oeis is None
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

    optimized.extend_to(2_000)
    reference.extend_to(2_000)

    assert optimized.terms == reference.terms
    assert optimized.used == set(optimized.terms)
    assert optimized.unused_multiplier_successors


@pytest.mark.parametrize("pair", [(2, 3), (2, 5), (3, 5)])
def test_every_return_after_two_free_terms_is_full(pair):
    p, q = pair
    generator = FullReturnEnotsWolleyGenerator(p=p, q=q)
    generator.extend_to(2_000)

    for index in range(2, len(generator.terms)):
        previous = generator.terms[index - 1]
        two_back = generator.terms[index - 2]
        value = generator.terms[index]
        if full_return_restriction_active(previous, two_back, p, q):
            has_p = value % p == 0
            has_q = value % q == 0
            assert has_p == has_q


def test_state_local_rejections_can_become_eligible_later():
    # These values are one-sided target candidates rejected in an earlier
    # two-free state, but they later occur once the restriction is inactive.
    cases = [
        ((2, 3), 56, 47),
        ((2, 5), 22, 19),
        ((3, 5), 6, 34),
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
