import pytest

from lex_earliest_seqs import registry
from lex_earliest_seqs.cache import open_run
from lex_earliest_seqs.zoo.factor_restricted_enots_wolley import (
    EWFactorPolicy,
    FactorRestrictedEnotsWolleyDefinition,
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
    "X000003": EWFactorPolicy(frozenset({3}), squarefree=False),
    "X000004": EWFactorPolicy(frozenset({2, 3}), squarefree=True),
    "X000005": EWFactorPolicy(frozenset({3}), squarefree=True),
}


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


@pytest.mark.parametrize(
    ("sequence_id", "alias", "prefix"),
    [
        (
            "X000001",
            "biprimary-ew",
            [
                1, 2, 6, 15, 35, 14, 12, 33, 55, 10,
                18, 21, 77, 22, 20, 45, 39, 26, 28, 63,
            ],
        ),
        (
            "X000002",
            "triprimary-ew",
            [
                1, 2, 6, 15, 35, 14, 12, 33, 55, 10,
                18, 21, 77, 22, 20, 45, 39, 26, 28, 63,
            ],
        ),
        (
            "X000003",
            "exact-triprimary-ew",
            [
                1, 2, 30, 105, 1001, 286, 60, 255, 1309, 154,
                78, 195, 385, 238, 102, 165, 455, 182, 66, 285,
            ],
        ),
        (
            "X000004",
            "squarefree-triprimary-ew",
            [
                1, 2, 6, 15, 35, 14, 22, 33, 21, 70,
                26, 39, 51, 34, 10, 55, 77, 42, 30, 65,
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
def test_factor_family_reference_prefixes(sequence_id, alias, prefix):
    definition = registry.resolve(sequence_id)
    assert registry.resolve(alias) is definition

    run = open_run(definition, use_cache=False)
    run.ensure(len(prefix))
    assert list(run.terms) == prefix


def test_multiplicity_flag_changes_biprimary_sequence_at_term_seven():
    nonsquarefree = open_run(registry.resolve("X000001"), use_cache=False)
    squarefree = open_run(registry.resolve("X000000"), use_cache=False)

    nonsquarefree.ensure(7)
    squarefree.ensure(7)

    assert nonsquarefree.terms[6] == 12
    assert squarefree.terms[6] == 22
    assert X000001_POLICY.allows(12)
    assert not X000000_POLICY.allows(12)


def test_squarefree_flag_is_independent_for_three_factor_supports():
    assert X000003_POLICY.allows(60)  # 2^2 * 3 * 5
    assert not X000005_POLICY.allows(60)
    assert X000003_POLICY.allows(30)
    assert X000005_POLICY.allows(30)

    assert X000002_POLICY.allows(12)  # two distinct prime factors
    assert not X000004_POLICY.allows(12)
    assert X000002_POLICY.allows(30)  # three distinct prime factors
    assert X000004_POLICY.allows(30)
