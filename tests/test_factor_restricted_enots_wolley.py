import pytest

from lex_earliest_seqs.zoo.factor_restricted_enots_wolley import (
    EWFactorPolicy,
    FactorRestrictedEnotsWolleyDefinition,
    FactorRestrictedEnotsWolleyGenerator,
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


def test_definition_factory_retains_policy_and_uses_generic_generator_by_default():
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
