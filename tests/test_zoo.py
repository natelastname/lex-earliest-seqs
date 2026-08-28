from lex_earliest_seqs import registry
from lex_earliest_seqs.cache import open_run


def test_enots_wolley_reference_prefix():
    definition = registry.resolve("ew")
    run = open_run(definition, use_cache=False)
    run.ensure(12)
    assert list(run.terms) == [1, 2, 6, 15, 35, 14, 12, 33, 55, 10, 18, 21]


def test_enots_wolley_prime_projection_is_registered():
    definition = registry.resolve("A336957")
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


def test_forced_squarefree_enots_wolley_reference_prefix():
    definition = registry.resolve("A000000")
    assert definition.oeis == "A000000"
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
