"""Binary Enots--Wolley sequence (OEIS A338833)."""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cache

from ..core import SequenceDefinition
from ..object_space import PositiveIntegers
from ..projections import binary_digit_projection


@cache
def bit_support(value: int) -> frozenset[int]:
    if value < 1:
        raise ValueError("value must be positive")
    return frozenset(bit for bit in range(value.bit_length()) if value & (1 << bit))


def is_candidate(value: int, previous: int, two_back: int) -> bool:
    support = bit_support(value)
    previous_support = bit_support(previous)
    two_back_support = bit_support(two_back)
    return (
        bool(support & previous_support)
        and not bool(support & two_back_support)
        and bool(support - previous_support)
    )


@dataclass
class BinaryEnotsWolleyGenerator:
    """Least-unused scan for the binary-support analogue of Enots--Wolley."""

    terms: list[int] = field(default_factory=lambda: [1, 3])
    used: set[int] = field(default_factory=lambda: {1, 3})
    smallest_unused: int = 2

    def extend_to(self, count: int) -> None:
        if count < 0:
            raise ValueError("count must be nonnegative")
        if count <= len(self.terms):
            return
        if len(self.terms) < 2:
            raise RuntimeError("BinaryEnotsWolleyGenerator state is missing initial terms")

        while len(self.terms) < count:
            candidate = self.smallest_unused
            while candidate in self.used or not is_candidate(
                candidate, self.terms[-1], self.terms[-2]
            ):
                candidate += 1
            self.terms.append(candidate)
            self.used.add(candidate)
            while self.smallest_unused in self.used:
                self.smallest_unused += 1


BINARY_ENOTS_WOLLEY = SequenceDefinition[int](
    id="A338833",
    oeis="A338833",
    name="Binary Enots--Wolley",
    aliases=("binary-ew", "bew"),
    generator_factory=BinaryEnotsWolleyGenerator,
    generator_version=1,
    definition_version=1,
    offset=1,
    object_space=PositiveIntegers(),
    projections={"binary-digits": binary_digit_projection()},
    description=(
        "Binary-support analogue of Enots--Wolley: consecutive terms share a set "
        "bit, terms two apart have disjoint set-bit supports, and each new term "
        "introduces a bit absent from its predecessor."
    ),
)
