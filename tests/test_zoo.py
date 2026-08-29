import pickle

from lex_earliest_seqs import registry
from lex_earliest_seqs.cache import open_run
from lex_earliest_seqs.zoo.binary_enots_wolley import (
    BinaryEnotsWolleyGenerator,
    bit_support,
    next_admissible_binary,
)
from lex_earliest_seqs.zoo.enots_wolley import EnotsWolleyGenerator, is_candidate
from lex_earliest_seqs.zoo.forced_squarefree_enots_wolley import (
    ForcedSquarefreeEnotsWolleyGenerator,
    is_squarefree,
)


def _reference_enots_wolley(count: int) -> list[int]:
    if count <= 2:
        return [1, 2][:count]

    terms = [1, 2]
    used = {1, 2}
    smallest_unused = 3
    while len(terms) < count:
        candidate = smallest_unused
        while candidate in used or not is_candidate(
            candidate, terms[-1], terms[-2]
        ):
            candidate += 1
        terms.append(candidate)
        used.add(candidate)
        while smallest_unused in used:
            smallest_unused += 1
    return terms


def _reference_binary_candidate(value: int, previous: int, two_back: int) -> bool:
    support = bit_support(value)
    previous_support = bit_support(previous)
    two_back_support = bit_support(two_back)
    return (
        bool(support & previous_support)
        and not bool(support & two_back_support)
        and bool(support - previous_support)
    )


def _reference_binary_enots_wolley(count: int) -> list[int]:
    if count <= 2:
        return [1, 3][:count]

    terms = [1, 3]
    used = {1, 3}
    smallest_unused = 2
    while len(terms) < count:
        candidate = smallest_unused
        while candidate in used or not _reference_binary_candidate(
            candidate, terms[-1], terms[-2]
        ):
            candidate += 1
        terms.append(candidate)
        used.add(candidate)
        while smallest_unused in used:
            smallest_unused += 1
    return terms


def _reference_forced_squarefree(count: int) -> list[int]:
    if count <= 2:
        return [1, 2][:count]

    terms = [1, 2]
    used = {1, 2}
    smallest_unused_squarefree = 3
    while len(terms) < count:
        candidate = smallest_unused_squarefree
        while (
            candidate in used
            or not is_squarefree(candidate)
            or not is_candidate(candidate, terms[-1], terms[-2])
        ):
            candidate += 1
        terms.append(candidate)
        used.add(candidate)
        while (
            smallest_unused_squarefree in used
            or not is_squarefree(smallest_unused_squarefree)
        ):
            smallest_unused_squarefree += 1
    return terms


def test_enots_wolley_reference_prefix():
    definition = registry.resolve("ew")
    run = open_run(definition, use_cache=False)
    run.ensure(12)
    assert list(run.terms) == [1, 2, 6, 15, 35, 14, 12, 33, 55, 10, 18, 21]


def test_enots_wolley_history_successors_match_direct_rule():
    generator = EnotsWolleyGenerator()
    generator.extend_to(2_000)
    assert generator.terms == _reference_enots_wolley(2_000)
    assert generator.used == set(generator.terms)
    assert generator.unused_multiplier_successors


def test_enots_wolley_successor_deletes_used_prime_multipliers():
    generator = EnotsWolleyGenerator()
    generator.used.update({4, 6, 8})

    assert generator._next_unused_multiplier(2, 1) == 5
    assert generator._find_multiplier_successor(2, 1) == 5
    assert generator.unused_multiplier_successors[2]


def test_enots_wolley_pickle_resumes_with_successors_without_radicals():
    generator = EnotsWolleyGenerator()
    generator.extend_to(500)
    assert generator.radicals is not None
    assert generator.unused_multiplier_successors
    successors = {
        prime: dict(parents)
        for prime, parents in generator.unused_multiplier_successors.items()
    }

    restored = pickle.loads(pickle.dumps(generator))
    assert restored.radicals is None
    assert restored.terms == generator.terms
    assert restored.used == generator.used
    assert restored.unused_multiplier_successors == successors
    assert restored.smallest_unused == generator.smallest_unused

    restored.extend_to(1_000)
    assert restored.terms == _reference_enots_wolley(1_000)


def test_enots_wolley_old_version3_pickle_lazily_adds_successor_maps():
    generator = EnotsWolleyGenerator()
    generator.extend_to(100)
    del generator.unused_multiplier_successors

    restored = pickle.loads(pickle.dumps(generator))
    assert not hasattr(restored, "unused_multiplier_successors")
    restored.extend_to(500)

    assert restored.unused_multiplier_successors
    assert restored.terms == _reference_enots_wolley(500)


def test_enots_wolley_prime_projection_is_registered():
    definition = registry.resolve("A336957")
    assert definition.generator_version == 3
    assert "prime-exponents" in definition.projections


def test_binary_successor_returns_exact_least_admissible_candidate():
    for previous, two_back, lower_bound in [
        (3, 1, 2),
        (6, 3, 4),
        (12, 6, 7),
        (17, 9, 2),
        (48, 20, 21),
    ]:
        expected = lower_bound
        while not _reference_binary_candidate(expected, previous, two_back):
            expected += 1
        assert next_admissible_binary(lower_bound, previous, two_back) == expected


def test_binary_enots_wolley_reference_prefix():
    definition = registry.resolve("A338833")
    run = open_run(definition, use_cache=False)
    run.ensure(15)
    assert list(run.terms) == [
        1,
        3,
        6,
        12,
        9,
        17,
        18,
        10,
        13,
        20,
        48,
        33,
        5,
        14,
        24,
    ]
    assert "binary-digits" in definition.projections


def test_binary_successor_generator_matches_direct_rule():
    generator = BinaryEnotsWolleyGenerator()
    generator.extend_to(1_000)
    assert generator.terms == _reference_binary_enots_wolley(1_000)


def test_binary_generator_pickle_resumes():
    generator = BinaryEnotsWolleyGenerator()
    generator.extend_to(500)
    restored = pickle.loads(pickle.dumps(generator))
    restored.extend_to(1_000)
    assert restored.terms == _reference_binary_enots_wolley(1_000)


def test_forced_squarefree_enots_wolley_reference_prefix():
    definition = registry.resolve("A399457")
    assert definition.oeis == "A399457"
    assert definition.generator_version == 3
    assert registry.resolve("forced-squarefree-ew") is definition

    run = open_run(definition, use_cache=False)
    run.ensure(20)
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
        70,
        26,
        39,
        51,
        34,
        10,
        55,
        77,
        42,
        30,
        65,
    ]
    assert "prime-exponents" in definition.projections


def test_forced_squarefree_history_frontiers_match_direct_rule():
    generator = ForcedSquarefreeEnotsWolleyGenerator()
    generator.extend_to(1_000)
    assert generator.terms == _reference_forced_squarefree(1_000)
    assert generator.cofactor_frontiers
    assert generator.cofactor_frontiers[2] > 2


def test_forced_squarefree_frontier_advances_past_used_products():
    generator = ForcedSquarefreeEnotsWolleyGenerator()
    generator.extend_to(500)

    before = generator._cofactor_frontier(2)
    generator.used.add(2 * before)
    after = generator._cofactor_frontier(2)

    assert after > before
    assert 2 * after not in generator.used


def test_forced_squarefree_pickle_resumes_with_frontiers_without_radicals():
    generator = ForcedSquarefreeEnotsWolleyGenerator()
    generator.extend_to(500)
    assert generator.radicals is not None
    assert generator.cofactor_frontiers
    frontiers = dict(generator.cofactor_frontiers)

    restored = pickle.loads(pickle.dumps(generator))
    assert restored.radicals is None
    assert restored.cofactor_frontiers == frontiers
    restored.extend_to(1_000)
    assert restored.terms == _reference_forced_squarefree(1_000)
