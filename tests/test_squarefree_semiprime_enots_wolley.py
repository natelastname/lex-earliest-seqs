from lex_earliest_seqs import registry
from lex_earliest_seqs.cache import open_run
from lex_earliest_seqs.zoo.enots_wolley import is_candidate
from lex_earliest_seqs.zoo.squarefree_semiprime_enots_wolley import (
    SquarefreeSemiprimeEnotsWolleyGenerator,
    is_squarefree_semiprime,
)


def _reference_squarefree_semiprime_enots_wolley(count: int) -> list[int]:
    if count <= 2:
        return [1, 2][:count]

    terms = [1, 2]
    used = {1, 2}
    while len(terms) < count:
        candidate = 6
        while (
            candidate in used
            or not is_squarefree_semiprime(candidate)
            or not is_candidate(candidate, terms[-1], terms[-2])
        ):
            candidate += 1
        terms.append(candidate)
        used.add(candidate)
    return terms


def test_squarefree_semiprime_enots_wolley_registry_and_prefix():
    definition = registry.resolve("X000000")
    assert definition.id == "squarefree-semiprime-enots-wolley"
    assert definition.oeis == "X000000"
    assert registry.resolve("semiprime-ew") is definition
    assert "prime-exponents" in definition.projections

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


def test_closed_form_generator_matches_direct_greedy_rule():
    generator = SquarefreeSemiprimeEnotsWolleyGenerator()
    generator.extend_to(80)
    assert generator.terms == _reference_squarefree_semiprime_enots_wolley(80)
