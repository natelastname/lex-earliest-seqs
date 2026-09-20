"""XOR-popcount sequence (OEIS A357578)."""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass, field

from ..core import SequenceDefinition
from ..object_space import PositiveIntegers
from ..projections import binary_digit_projection


def next_same_hamming_weight(value: int) -> int:
    """Return the least integer greater than ``value`` with the same Hamming weight.

    The update rolls the lowest block of 1-bits one place to the left, then pulls
    the displaced 1-bits back down as far as possible. Repeated calls therefore
    walk a Hamming-weight class in ordinary increasing order without scanning
    integers from other weight classes.
    """

    if value < 1:
        raise ValueError("value must be positive")

    # Roll over the first (lowest) block of ones.
    lowest_one = value & -value
    rolled = value + lowest_one

    # Pull the displaced leading ones back down behind the rolled bit.
    pulled_back = value >> lowest_one.bit_length()
    pulled_back &= ~(pulled_back + 1)

    return rolled | pulled_back


def hamming_weight_class(weight: int) -> Iterator[int]:
    """Yield positive integers of Hamming weight ``weight`` in increasing order."""

    if weight < 1:
        raise ValueError("weight must be positive")

    value = (1 << weight) - 1
    while True:
        yield value
        value = next_same_hamming_weight(value)


@dataclass
class XorpopGenerator:
    """FIFO generator for A357578 using one frontier per Hamming-weight class."""

    terms: list[int] = field(default_factory=lambda: [1, 2])
    weight_frontiers: dict[int, int] = field(default_factory=lambda: {1: 4})

    def extend_to(self, count: int) -> None:
        if count < 0:
            raise ValueError("count must be nonnegative")
        if count <= len(self.terms):
            return
        if len(self.terms) < 2:
            raise RuntimeError("XorpopGenerator state is missing initial terms")

        while len(self.terms) < count:
            weight = (self.terms[-1] ^ self.terms[-2]).bit_count()
            candidate = self.weight_frontiers.get(weight, (1 << weight) - 1)

            self.terms.append(candidate)
            self.weight_frontiers[weight] = next_same_hamming_weight(candidate)


XORPOP = SequenceDefinition[int](
    id="A357578",
    oeis="A357578",
    name="XOR-popcount",
    aliases=("xorpop", "xor-pop"),
    generator_factory=XorpopGenerator,
    generator_version=1,
    definition_version=1,
    offset=1,
    object_space=PositiveIntegers(),
    projections={"binary-digits": binary_digit_projection()},
    description=(
        "Lexicographically earliest sequence of distinct positive integers in which "
        "the Hamming weight of each new term equals the Hamming weight of the XOR "
        "of the previous two terms. The generator uses the FIFO order within each "
        "Hamming-weight class."
    ),
)
