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

    This is a digit-DP successor rather than an integer scan. A valid candidate
    must avoid every set bit of ``two_back``, use at least one set bit of
    ``previous`` that is not forbidden by ``two_back``, and use at least one bit
    absent from ``previous``. The DP has only the bit position and three Boolean
    flags as state, so finding the successor takes O(log candidate) work.
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

    # One bit beyond all inputs is always available as a new, un-forbidden bit,
    # so this fixed width is enough whenever an admissible successor exists.
    width = max(
        lower_bound.bit_length(),
        previous.bit_length(),
        two_back.bit_length(),
    ) + 1

    @cache
    def feasible(
        position: int,
        already_greater: bool,
        has_shared_bit: bool,
        has_new_bit: bool,
    ) -> bool:
        if position < 0:
            return has_shared_bit and has_new_bit

        lower_bit = (lower_bound >> position) & 1
        forbidden = (two_back >> position) & 1
        previous_bit = (previous >> position) & 1
        shared_bit = (shared_mask >> position) & 1

        for bit in (0, 1):
            if forbidden and bit:
                continue
            if not already_greater and bit < lower_bit:
                continue

            next_greater = already_greater or bit > lower_bit
            next_shared = has_shared_bit or bool(bit and shared_bit)
            next_new = has_new_bit or bool(bit and not previous_bit)
            if feasible(position - 1, next_greater, next_shared, next_new):
                return True
        return False

    if not feasible(width - 1, False, False, False):
        raise RuntimeError("binary EW successor DP found no admissible candidate")

    result = 0
    already_greater = False
    has_shared_bit = False
    has_new_bit = False

    for position in range(width - 1, -1, -1):
        lower_bit = (lower_bound >> position) & 1
        forbidden = (two_back >> position) & 1
        previous_bit = (previous >> position) & 1
        shared_bit = (shared_mask >> position) & 1

        for bit in (0, 1):
            if forbidden and bit:
                continue
            if not already_greater and bit < lower_bit:
                continue

            next_greater = already_greater or bit > lower_bit
            next_shared = has_shared_bit or bool(bit and shared_bit)
            next_new = has_new_bit or bool(bit and not previous_bit)
            if not feasible(position - 1, next_greater, next_shared, next_new):
                continue

            if bit:
                result |= 1 << position
            already_greater = next_greater
            has_shared_bit = next_shared
            has_new_bit = next_new
            break
        else:  # pragma: no cover - guarded by the feasibility pass above
            raise RuntimeError("failed to reconstruct binary EW successor")

    return result


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
