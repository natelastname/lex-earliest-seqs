import pickle

from lex_earliest_seqs import registry
from lex_earliest_seqs.cache import open_run
from lex_earliest_seqs.zoo.binary_enots_wolley import (
    BinaryEnotsWolleyGenerator,
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


def _reference_bit_support(value: int) -> set[int]:
    return {bit for bit in range(value.bit_length()) if value & (1 << bit)}


def _reference_binary_candidate(value: int, previous: int, two_back: int) -> bool:
    support = _reference_bit_support(value)
    previous_support = _reference_bit_support(previous)
    two_back_support = _reference_bit_support(two_back)
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


def _reference_forced_squarefree_enots_wolley(count: int) -> list[int]:
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


def test_fast_enots_wolley_matches_direct_rule():
    generator = EnotsWolleyGenerator()
    generator.extend_to(1_000)
    assert generator.terms == _reference_enots_wolley(1_000)
    assert generator.used == set(generator.terms)


def test_fast_enots_wolley_pickle_resumes_without_persisting_radicals():
    generator = EnotsWolleyGenerator()
    generator.extend_to(500)
    assert generator.radicals is not None

    restored = pickle.loads(pickle.dumps(generator))
    assert restored.radicals is None
    assert restored.terms == generator.terms
    assert restored.used == generator.used
    assert isinstance(restored.used, set)
    assert restored.smallest_unused == generator.smallest_unused

    restored.extend_to(1_000)
    assert restored.terms == _reference_enots_wolley(1_000)


def test_enots_wolley_prime_projection_is_registered():
    definition = registry.resolve("A336957")
    assert definition.generator_version == 3
    assert "prime-exponents" in definition.projections


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


def test_binary_successor_is_least_admissible_at_generated_states():
    terms = _reference_binary_enots_wolley(200)
    used: set[int] = set()
    smallest_unused = 1

    for index, value in enumerate(terms):
        used.add(value)
        while smallest_unused in used:
            smallest_unused += 1
        if index < 1 or index + 1 >= len(terms):
            continue

        successor = next_admissible_binary(
            smallest_unused,
            terms[index],
            terms[index - 1],
        )
        brute_force = smallest_unused
        while not _reference_binary_candidate(
            brute_force,
            terms[index],
            terms[index - 1],
        ):
            brute_force += 1
        assert successor == brute_force


def test_forced_squarefree_enots_wolley_reference_prefix():
    definition = registry.resolve("A399457")
    assert definition.oeis == "A399457"
    assert definition.generator_version == 2
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


def test_squarefree_stream_generator_matches_direct_rule():
    generator = ForcedSquarefreeEnotsWolleyGenerator()
    generator.extend_to(500)
    assert generator.terms == _reference_forced_squarefree_enots_wolley(500)


def test_squarefree_stream_generator_pickle_resumes_without_radicals():
    expected = _reference_forced_squarefree_enots_wolley(500)
    generator = ForcedSquarefreeEnotsWolleyGenerator()
    generator.extend_to(300)
    assert generator.radicals is not None

    restored = pickle.loads(pickle.dumps(generator))
    assert restored.radicals is None
    assert restored.terms == generator.terms
    assert restored.used == generator.used
    assert restored.smallest_unused_squarefree == generator.smallest_unused_squarefree

    restored.extend_to(500)
    assert restored.terms == expected
