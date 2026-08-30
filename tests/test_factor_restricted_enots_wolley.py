import pytest

from lex_earliest_seqs.zoo.factor_restricted_enots_wolley import (
    EWFactorPolicy,
    FactorRestrictedEnotsWolleyDefinition,
    FactorRestrictedEnotsWolleyGenerator,
    big_omega,
    make_factor_restricted_enots_wolley_definition,
)


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


def test_squarefree_is_independent_of_allowed_omega():
    multiplicity_allowed = EWFactorPolicy(
        allowed_omega=frozenset({3}),
        squarefree=False,
    )
    squarefree_only = EWFactorPolicy(
        allowed_omega=frozenset({3}),
        squarefree=True,
    )

    assert multiplicity_allowed.allows(12)
    assert multiplicity_allowed.allows(30)
    assert squarefree_only.allows(30)
    assert not squarefree_only.allows(12)


def test_definition_factory_retains_policy_and_uses_generic_generator_by_default():
    policy = EWFactorPolicy(
        allowed_omega=frozenset({2, 3}),
        squarefree=False,
    )
    definition = make_factor_restricted_enots_wolley_definition(
        id="X999999",
        name="Test finite-Omega EW",
        policy=policy,
    )

    assert isinstance(definition, FactorRestrictedEnotsWolleyDefinition)
    assert definition.factor_policy == policy
    generator = definition.generator_factory()
    assert isinstance(generator, FactorRestrictedEnotsWolleyGenerator)
    assert generator.policy == policy
    assert "prime-exponents" in definition.projections


def test_omega_two_ew_is_unchanged_by_squarefree_flag():
    multiplicity_allowed = FactorRestrictedEnotsWolleyGenerator(
        policy=EWFactorPolicy(frozenset({2}), squarefree=False)
    )
    squarefree_only = FactorRestrictedEnotsWolleyGenerator(
        policy=EWFactorPolicy(frozenset({2}), squarefree=True)
    )

    multiplicity_allowed.extend_to(80)
    squarefree_only.extend_to(80)

    assert multiplicity_allowed.terms == squarefree_only.terms
