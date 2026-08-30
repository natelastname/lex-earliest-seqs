import pytest

from lex_earliest_seqs.projections import prime_factorization
from lex_earliest_seqs.zoo.factor_restricted_enots_wolley import (
    EWFactorPolicy,
    FactorRestrictedEnotsWolleyDefinition,
    FactorRestrictedEnotsWolleyGenerator,
    ReferenceFactorRestrictedEnotsWolleyGenerator,
    big_omega,
    make_factor_restricted_enots_wolley_definition,
    omega,
)


def test_omega_counts_distinct_prime_factors():
    assert omega(1) == 0
    assert omega(12) == 2
    assert omega(18) == 2
    assert omega(30) == 3


def test_big_omega_counts_prime_factors_with_multiplicity():
    assert big_omega(1) == 0
    assert big_omega(12) == 3
    assert big_omega(18) == 3
    assert big_omega(30) == 3


def test_factor_policy_normalizes_and_validates_allowed_omega():
    policy = EWFactorPolicy(allowed_omega={2, 3})
    assert policy.allowed_omega == frozenset({2, 3})

    with pytest.raises(ValueError):
        EWFactorPolicy(allowed_omega=frozenset())
    with pytest.raises(ValueError):
        EWFactorPolicy(allowed_omega=frozenset({-1, 2}))
    with pytest.raises(TypeError):
        EWFactorPolicy(allowed_omega=frozenset({2, 3.0}))


def test_squarefree_is_independent_of_distinct_prime_count():
    multiplicity_allowed = EWFactorPolicy(
        allowed_omega=frozenset({2}),
        squarefree=False,
    )
    squarefree_only = EWFactorPolicy(
        allowed_omega=frozenset({2}),
        squarefree=True,
    )

    assert multiplicity_allowed.allows(12)
    assert multiplicity_allowed.allows(18)
    assert multiplicity_allowed.allows(6)
    assert not squarefree_only.allows(12)
    assert not squarefree_only.allows(18)
    assert squarefree_only.allows(6)


def test_definition_factory_retains_policy_and_uses_optimized_generator_by_default():
    policy = EWFactorPolicy(
        allowed_omega=frozenset({2, 3}),
        squarefree=False,
    )
    definition = make_factor_restricted_enots_wolley_definition(
        id="X999999",
        name="Test finite-omega EW",
        policy=policy,
    )

    assert isinstance(definition, FactorRestrictedEnotsWolleyDefinition)
    assert definition.factor_policy == policy
    generator = definition.generator_factory()
    assert isinstance(generator, FactorRestrictedEnotsWolleyGenerator)
    assert generator.policy == policy
    assert "prime-exponents" in definition.projections


def test_reference_generator_remains_available_as_independent_oracle():
    policy = EWFactorPolicy(frozenset({2}), squarefree=False)
    reference = ReferenceFactorRestrictedEnotsWolleyGenerator(policy=policy)
    reference.extend_to(10)
    assert reference.terms == [1, 2, 6, 15, 35, 14, 12, 33, 55, 10]


def test_biprimary_multiplicity_allowed_differs_from_squarefree_variant():
    multiplicity_allowed = FactorRestrictedEnotsWolleyGenerator(
        policy=EWFactorPolicy(frozenset({2}), squarefree=False)
    )
    squarefree_only = FactorRestrictedEnotsWolleyGenerator(
        policy=EWFactorPolicy(frozenset({2}), squarefree=True)
    )

    multiplicity_allowed.extend_to(7)
    squarefree_only.extend_to(7)

    assert multiplicity_allowed.terms == [1, 2, 6, 15, 35, 14, 12]
    assert squarefree_only.terms == [1, 2, 6, 15, 35, 14, 22]


def test_persistent_successor_map_deletes_policy_failures_and_used_products():
    generator = FactorRestrictedEnotsWolleyGenerator(
        policy=EWFactorPolicy(frozenset({2}), squarefree=True)
    )

    # For stream prime 2, multipliers 1 and 2 produce 2 and 4, which have
    # omega=1. The first globally policy-eligible product is 2*3=6.
    assert generator._next_persistently_eligible_multiplier(2, 1) == 3
    parents = generator.multiplier_successors[2]
    assert 1 in parents
    assert 2 in parents

    # Once 6 is globally used, multiplier 3 is permanently deletable as well;
    # multiplier 4 produces nonsquarefree 8, so the next head is 2*5=10.
    generator.used.add(6)
    assert generator._next_persistently_eligible_multiplier(2, 1) == 5
    assert 3 in parents
    assert 4 in parents


def test_successor_find_path_compresses_deleted_multiplier_chain():
    generator = FactorRestrictedEnotsWolleyGenerator(
        policy=EWFactorPolicy(frozenset({2}), squarefree=False)
    )
    generator.multiplier_successors[2] = {1: 2, 2: 3, 3: 7}

    assert generator._find_multiplier_successor(2, 1) == 7
    assert generator.multiplier_successors[2][1] == 7
    assert generator.multiplier_successors[2][2] == 7
    assert generator.multiplier_successors[2][3] == 7


def test_delete_multiplier_is_idempotent_and_links_to_live_successor():
    generator = FactorRestrictedEnotsWolleyGenerator(
        policy=EWFactorPolicy(frozenset({2, 3}), squarefree=False)
    )

    assert generator._delete_multiplier(3, 20)
    assert generator._find_multiplier_successor(3, 20) == 21
    assert not generator._delete_multiplier(3, 20)

    generator._delete_multiplier(3, 21)
    assert generator._find_multiplier_successor(3, 20) == 22
    assert generator.multiplier_successors[3][20] == 22


def test_retire_used_value_deletes_every_prime_stream_representation():
    generator = FactorRestrictedEnotsWolleyGenerator(
        policy=EWFactorPolicy(frozenset({2, 3}), squarefree=False)
    )
    generator.used.add(60)

    # 60 = 2*30 = 3*20 = 5*12. All three stream representations are
    # permanently dead as soon as 60 has been selected.
    assert generator._retire_used_value(60) == 3
    for stream_prime, multiplier in ((2, 30), (3, 20), (5, 12)):
        assert generator._find_multiplier_successor(stream_prime, multiplier) > multiplier

    # Re-retirement is a no-op rather than creating a second deletion chain.
    assert generator._retire_used_value(60) == 0


def test_extend_eagerly_retires_all_representations_of_selected_terms():
    generator = FactorRestrictedEnotsWolleyGenerator(
        policy=EWFactorPolicy(frozenset({2, 3}), squarefree=False)
    )
    generator.extend_to(100)

    for value in generator.terms[2:]:
        for stream_prime, _ in prime_factorization(value):
            multiplier = value // stream_prime
            assert (
                generator._find_multiplier_successor(stream_prime, multiplier)
                > multiplier
            )


def test_eager_retirement_does_not_change_sequence_against_reference():
    policy = EWFactorPolicy(frozenset({2, 3}), squarefree=False)
    optimized = FactorRestrictedEnotsWolleyGenerator(policy=policy)
    reference = ReferenceFactorRestrictedEnotsWolleyGenerator(policy=policy)

    optimized.extend_to(1_000)
    reference.extend_to(1_000)

    assert optimized.terms == reference.terms
