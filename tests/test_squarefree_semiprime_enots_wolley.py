from lex_earliest_seqs import registry
from lex_earliest_seqs.cache import open_run
from lex_earliest_seqs.zoo.factor_restricted_enots_wolley import (
    FactorRestrictedEnotsWolleyDefinition,
    FactorRestrictedEnotsWolleyGenerator,
)
from lex_earliest_seqs.zoo.squarefree_semiprime_enots_wolley import (
    X000000_POLICY,
    SquarefreeSemiprimeEnotsWolleyGenerator,
    is_squarefree_semiprime,
)


def test_squarefree_semiprime_enots_wolley_registry_and_prefix():
    definition = registry.resolve("X000000")
    assert isinstance(definition, FactorRestrictedEnotsWolleyDefinition)
    assert definition.id == "X000000"
    assert definition.oeis is None
    assert definition.generator_version == 2
    assert definition.factor_policy == X000000_POLICY
    assert registry.resolve("semiprime-ew") is definition
    assert registry.resolve("squarefree-semiprime-enots-wolley") is definition
    assert "prime-exponents" in definition.projections

    generator = definition.generator_factory()
    assert isinstance(generator, SquarefreeSemiprimeEnotsWolleyGenerator)
    assert generator.policy == X000000_POLICY
    assert X000000_POLICY.allowed_omega == frozenset({2})
    assert X000000_POLICY.squarefree

    run = open_run(definition, use_cache=False)
    run.ensure(35)
    assert list(run.terms) == [
        1,
        2,
        6,
        15,
        35,
        14,
        22,
        33,
        21,
        91,
        26,
        10,
        55,
        77,
        119,
        34,
        38,
        57,
        39,
        65,
        85,
        51,
        69,
        46,
        58,
        87,
        93,
        62,
        74,
        111,
        123,
        82,
        86,
        129,
        141,
    ]


def test_x000000_closed_form_matches_generic_factor_restricted_generator():
    optimized = SquarefreeSemiprimeEnotsWolleyGenerator()
    reference = FactorRestrictedEnotsWolleyGenerator(policy=X000000_POLICY)

    optimized.extend_to(80)
    reference.extend_to(80)

    assert optimized.terms == reference.terms


def test_squarefree_semiprime_predicate_is_policy_backed():
    assert is_squarefree_semiprime(6)
    assert is_squarefree_semiprime(35)
    assert not is_squarefree_semiprime(4)
    assert not is_squarefree_semiprime(12)
