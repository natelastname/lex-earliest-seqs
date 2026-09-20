import pickle
from itertools import islice

from lex_earliest_seqs import registry
from lex_earliest_seqs.cache import open_run
from lex_earliest_seqs.zoo.xorpop import (
    XorpopGenerator,
    hamming_weight_class,
    next_same_hamming_weight,
)


OEIS_PREFIX = [
    1,
    2,
    3,
    4,
    7,
    5,
    8,
    11,
    6,
    13,
    14,
    9,
    19,
    21,
    10,
    31,
    22,
    12,
    25,
    26,
    17,
    28,
    35,
    63,
    37,
    38,
    18,
    41,
    47,
    20,
    55,
    42,
    15,
    44,
    49,
    23,
    50,
    52,
    24,
    56,
    16,
    33,
    67,
    69,
    34,
    59,
    70,
    95,
    73,
    74,
    36,
    61,
    76,
    27,
    62,
    81,
    111,
    79,
    32,
    119,
    87,
    64,
    29,
    91,
    82,
    40,
    93,
    94,
    48,
    103,
    107,
    65,
    84,
    88,
]


def test_same_weight_successor_walks_weight_class_in_natural_order():
    assert list(islice(hamming_weight_class(2), 10)) == [
        3,
        5,
        6,
        9,
        10,
        12,
        17,
        18,
        20,
        24,
    ]

    value = 7
    successors = []
    for _ in range(9):
        successors.append(value)
        value = next_same_hamming_weight(value)

    assert successors == [7, 11, 13, 14, 19, 21, 22, 25, 26]
    assert all(item.bit_count() == 3 for item in successors)


def test_xorpop_matches_oeis_prefix_and_registry_metadata():
    definition = registry.resolve("A357578")
    assert definition.oeis == "A357578"
    assert registry.resolve("xorpop") is definition
    assert "binary-digits" in definition.projections

    run = open_run(definition, use_cache=False)
    run.ensure(len(OEIS_PREFIX))

    assert list(run.terms) == OEIS_PREFIX


def test_xorpop_fifo_law_and_recurrence_on_long_prefix():
    generator = XorpopGenerator()
    generator.extend_to(10_000)

    assert len(set(generator.terms)) == len(generator.terms)

    for index in range(2, len(generator.terms)):
        assert generator.terms[index].bit_count() == (
            generator.terms[index - 1] ^ generator.terms[index - 2]
        ).bit_count()

    terms_by_weight: dict[int, list[int]] = {}
    for term in generator.terms:
        terms_by_weight.setdefault(term.bit_count(), []).append(term)

    for weight, terms in terms_by_weight.items():
        assert terms == list(islice(hamming_weight_class(weight), len(terms)))
        assert generator.weight_frontiers[weight] == next_same_hamming_weight(terms[-1])


def test_xorpop_generator_pickle_resumes():
    generator = XorpopGenerator()
    generator.extend_to(1_000)

    restored = pickle.loads(pickle.dumps(generator))
    restored.extend_to(2_000)

    expected = XorpopGenerator()
    expected.extend_to(2_000)

    assert restored.terms == expected.terms
    assert restored.weight_frontiers == expected.weight_frontiers
