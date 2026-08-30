import pickle
from collections import defaultdict

import pytest

from lex_earliest_seqs import registry
from lex_earliest_seqs.cache import open_run
from lex_earliest_seqs.zoo.enots_wolley import prime_support
from lex_earliest_seqs.zoo.factor_restricted_enots_wolley import (
    EWFactorPolicy,
    FactorRestrictedEnotsWolleyDefinition,
    FactorRestrictedEnotsWolleyGenerator,
    ReferenceFactorRestrictedEnotsWolleyGenerator,
)
from lex_earliest_seqs.zoo.primary_enots_wolley import (
    X000001_POLICY,
    X000002_POLICY,
    X000003_POLICY,
    X000004_POLICY,
    X000005_POLICY,
)
from lex_earliest_seqs.zoo.squarefree_semiprime_enots_wolley import X000000_POLICY


EXPECTED_POLICIES = {
    "X000000": EWFactorPolicy(frozenset({2}), squarefree=True),
    "X000001": EWFactorPolicy(frozenset({2}), squarefree=False),
    "X000002": EWFactorPolicy(frozenset({2, 3}), squarefree=False),
    "X000003": EWFactorPolicy(frozenset({2, 3}), squarefree=True),
    "X000004": EWFactorPolicy(frozenset({3}), squarefree=False),
    "X000005": EWFactorPolicy(frozenset({3}), squarefree=True),
}

GENERIC_IDS = ("X000001", "X000002", "X000003", "X000004", "X000005")


@pytest.mark.parametrize(
    ("sequence_id", "policy"),
    [
        ("X000000", X000000_POLICY),
        ("X000001", X000001_POLICY),
        ("X000002", X000002_POLICY),
        ("X000003", X000003_POLICY),
        ("X000004", X000004_POLICY),
        ("X000005", X000005_POLICY),
    ],
)
def test_all_six_factor_family_members_are_registered(sequence_id, policy):
    definition = registry.resolve(sequence_id)
    assert isinstance(definition, FactorRestrictedEnotsWolleyDefinition)
    assert definition.oeis is None
    assert definition.factor_policy == policy == EXPECTED_POLICIES[sequence_id]
    assert "prime-exponents" in definition.projections
    if sequence_id != "X000000":
        assert definition.generator_version == 2
        assert isinstance(definition.generator_factory(), FactorRestrictedEnotsWolleyGenerator)


@pytest.mark.parametrize(
    ("sequence_id", "alias", "prefix"),
    [
        (
            "X000001",
            "biprimary-ew",
            [
                1, 2, 6, 15, 35, 14, 12, 33, 55, 10,
                18, 21, 77, 22, 20, 45, 39, 26, 28, 63,
                51, 34, 38, 57, 69, 46, 40, 65, 91, 56,
            ],
        ),
        (
            "X000002",
            "triprimary-ew",
            [
                1, 2, 6, 15, 35, 14, 12, 33, 55, 10,
                18, 21, 77, 22, 20, 45, 39, 26, 28, 63,
                51, 34, 38, 57, 69, 46, 40, 65, 91, 42,
            ],
        ),
        (
            "X000003",
            "squarefree-triprimary-ew",
            [
                1, 2, 6, 15, 35, 14, 22, 33, 21, 70,
                26, 39, 51, 34, 10, 55, 77, 42, 30, 65,
            ],
        ),
        (
            "X000004",
            "exact-triprimary-ew",
            [
                1, 2, 30, 105, 1001, 286, 60, 255, 1309, 154,
                78, 195, 385, 238, 102, 165, 455, 182, 66, 285,
            ],
        ),
        (
            "X000005",
            "squarefree-exact-triprimary-ew",
            [
                1, 2, 30, 105, 1001, 286, 102, 255, 385, 154,
                78, 195, 595, 238, 66, 165, 455, 182, 114, 285,
            ],
        ),
    ],
)
def test_factor_family_regression_prefixes(sequence_id, alias, prefix):
    definition = registry.resolve(sequence_id)
    assert registry.resolve(alias) is definition

    run = open_run(definition, use_cache=False)
    run.ensure(len(prefix))
    assert list(run.terms) == prefix


@pytest.mark.parametrize("sequence_id", GENERIC_IDS)
def test_optimized_generator_matches_independent_reference_through_1000_terms(
    sequence_id,
):
    definition = registry.resolve(sequence_id)
    policy = definition.factor_policy
    assert policy is not None

    optimized = FactorRestrictedEnotsWolleyGenerator(policy=policy)
    reference = ReferenceFactorRestrictedEnotsWolleyGenerator(policy=policy)
    optimized.extend_to(1_000)
    reference.extend_to(1_000)

    assert optimized.terms == reference.terms
    assert optimized.used == reference.used
    assert optimized.multiplier_successors


@pytest.mark.parametrize("sequence_id", GENERIC_IDS)
def test_optimized_generator_pickle_resume_matches_reference(sequence_id):
    definition = registry.resolve(sequence_id)
    policy = definition.factor_policy
    assert policy is not None

    optimized = FactorRestrictedEnotsWolleyGenerator(policy=policy)
    optimized.extend_to(250)
    successor_snapshot = {
        prime: dict(parents)
        for prime, parents in optimized.multiplier_successors.items()
    }

    restored = pickle.loads(pickle.dumps(optimized))
    assert restored.terms == optimized.terms
    assert restored.used == optimized.used
    assert restored.multiplier_successors == successor_snapshot

    restored.extend_to(750)
    reference = ReferenceFactorRestrictedEnotsWolleyGenerator(policy=policy)
    reference.extend_to(750)
    assert restored.terms == reference.terms


@pytest.mark.parametrize("sequence_id", GENERIC_IDS)
def test_generated_prefix_satisfies_full_definition_and_is_injective(sequence_id):
    definition = registry.resolve(sequence_id)
    policy = definition.factor_policy
    assert policy is not None

    generator = FactorRestrictedEnotsWolleyGenerator(policy=policy)
    generator.extend_to(1_000)
    assert len(generator.terms) == len(set(generator.terms))

    for index in range(2, len(generator.terms)):
        value = generator.terms[index]
        previous = generator.terms[index - 1]
        two_back = generator.terms[index - 2]
        support = prime_support(value)
        previous_support = prime_support(previous)

        assert policy.allows(value)
        assert support & previous_support
        assert not support & prime_support(two_back)
        assert support - previous_support


@pytest.mark.parametrize("sequence_id", ("X000001", "X000002", "X000004"))
def test_each_reused_exact_support_is_served_in_numerical_queue_order(sequence_id):
    """Selected values on one support must form an initial exact-support queue."""

    definition = registry.resolve(sequence_id)
    generator = definition.generator_factory()
    generator.extend_to(500)

    by_support = defaultdict(list)
    for value in generator.terms[2:]:
        by_support[prime_support(value)].append(value)

    checked = 0
    for support, values in by_support.items():
        if len(values) < 2:
            continue
        last = values[-1]
        exact_support_prefix = [
            value
            for value in range(2, last + 1)
            if prime_support(value) == support
        ]
        assert values == exact_support_prefix
        checked += 1
        if checked >= 20:
            break

    assert checked > 0


@pytest.mark.parametrize("sequence_id", ("X000003", "X000005"))
def test_squarefree_variants_never_reuse_an_exact_prime_support(sequence_id):
    definition = registry.resolve(sequence_id)
    generator = definition.generator_factory()
    generator.extend_to(500)

    supports = [prime_support(value) for value in generator.terms[2:]]
    assert len(supports) == len(set(supports))


def test_multiplicity_flag_changes_biprimary_sequence_at_term_seven():
    nonsquarefree = open_run(registry.resolve("X000001"), use_cache=False)
    squarefree = open_run(registry.resolve("X000000"), use_cache=False)

    nonsquarefree.ensure(7)
    squarefree.ensure(7)

    assert nonsquarefree.terms[6] == 12
    assert squarefree.terms[6] == 22
    assert X000001_POLICY.allows(12)
    assert not X000000_POLICY.allows(12)


def test_allowing_three_factor_supports_first_changes_nonsquarefree_term_30():
    biprimary = open_run(registry.resolve("X000001"), use_cache=False)
    triprimary = open_run(registry.resolve("X000002"), use_cache=False)

    biprimary.ensure(30)
    triprimary.ensure(30)

    assert list(biprimary.terms[:29]) == list(triprimary.terms[:29])
    assert biprimary.terms[29] == 56
    assert triprimary.terms[29] == 42


def test_squarefree_flag_is_independent_for_three_factor_supports():
    assert X000004_POLICY.allows(60)  # 2^2 * 3 * 5
    assert not X000005_POLICY.allows(60)
    assert X000004_POLICY.allows(30)
    assert X000005_POLICY.allows(30)

    assert X000002_POLICY.allows(12)  # two distinct prime factors
    assert not X000003_POLICY.allows(12)
    assert X000002_POLICY.allows(30)  # three distinct prime factors
    assert X000003_POLICY.allows(30)
