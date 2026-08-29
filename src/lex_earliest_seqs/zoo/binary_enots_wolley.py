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
    """Return whether ``value`` satisfies the binary EW adjacency rule."""

    if value < 1:
        return False
    return (
        (value & previous) != 0
        and (value & two_back) == 0
        and (value & ~previous) != 0
    )


def next_admissible_binary(
    lower_bound: int,
    previous: int,
    two_back: int,
) -> int:
    """Return the least binary-EW candidate at least ``lower_bound``.

    No intervening integers are tested. If ``lower_bound`` itself is invalid,
    the successor must have a highest bit where it differs from the lower bound;
    at that pivot the lower-bound bit is 0 and the successor bit is 1. We scan
    possible pivots from least to most significant. For a chosen pivot, all more
    significant bits remain equal to the lower bound, while the less significant
    suffix is filled minimally with whichever of the two required bit classes is
    still missing:

    * a bit shared with ``previous`` but absent from ``two_back``;
    * a bit absent from both ``previous`` and ``two_back``.

    The first feasible pivot therefore gives the least admissible successor in
    O(log candidate) bit operations.
    """

    if lower_bound < 1:
        raise ValueError("lower_bound must be positive")
    if previous < 1 or two_back < 1:
        raise ValueError("previous and two_back must be positive")

    shared_mask = previous & ~two_back
    if shared_mask == 0:
        raise RuntimeError(
            "binary EW state has no bit available to share with the predecessor"
        )

    if is_candidate(lower_bound, previous, two_back):
        return lower_bound

    width = max(
        lower_bound.bit_length(),
        previous.bit_length(),
        two_back.bit_length(),
    ) + 1

    for pivot in range(width):
        # The highest differing bit of a larger number must change 0 -> 1, and
        # that bit may not be forbidden by the two-back support.
        if (lower_bound >> pivot) & 1:
            continue
        if (two_back >> pivot) & 1:
            continue

        # Bits above the pivot stay exactly equal to lower_bound. Any forbidden
        # set bit there makes this pivot impossible.
        upper = lower_bound >> (pivot + 1)
        if upper & (two_back >> (pivot + 1)):
            continue

        has_shared = bool(upper & (shared_mask >> (pivot + 1))) or bool(
            (shared_mask >> pivot) & 1
        )
        has_new = bool(upper & ~(previous >> (pivot + 1))) or not bool(
            (previous >> pivot) & 1
        )

        lower_mask = (1 << pivot) - 1
        suffix = 0

        if not has_shared:
            available_shared = shared_mask & lower_mask
            if not available_shared:
                continue
            suffix |= available_shared & -available_shared

        if not has_new:
            available_new = ~(previous | two_back) & lower_mask
            if not available_new:
                continue
            suffix |= available_new & -available_new

        return (
            (lower_bound >> (pivot + 1) << (pivot + 1))
            | (1 << pivot)
            | suffix
        )

    # The extra top bit included in ``width`` is always a new, un-forbidden bit;
    # with ``shared_mask != 0`` a feasible pivot must therefore exist.
    raise RuntimeError("binary EW successor construction found no candidate")


@dataclass
class BinaryEnotsWolleyGenerator:
    """Binary EW generator using an exact admissible-candidate successor."""

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
            candidate = next_admissible_binary(
                self.smallest_unused,
                self.terms[-1],
                self.terms[-2],
            )
            while candidate in self.used:
                candidate = next_admissible_binary(
                    candidate + 1,
                    self.terms[-1],
                    self.terms[-2],
                )

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
