"""Enots--Wolley sequence definition (OEIS A336957)."""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cache
from math import gcd

from ..core import SequenceDefinition
from ..object_space import PositiveIntegers
from ..projections import prime_exponent_projection


@cache
def prime_support(value: int) -> frozenset[int]:
    if value < 1:
        raise ValueError("value must be positive")
    factors: set[int] = set()
    remaining = value
    divisor = 2
    while divisor * divisor <= remaining:
        if remaining % divisor == 0:
            factors.add(divisor)
            while remaining % divisor == 0:
                remaining //= divisor
        divisor = 3 if divisor == 2 else divisor + 2
    if remaining > 1:
        factors.add(remaining)
    return frozenset(factors)


def is_candidate(value: int, previous: int, two_back: int) -> bool:
    return (
        value >= 1
        and gcd(value, previous) != 1
        and gcd(value, two_back) == 1
        and bool(prime_support(value) - prime_support(previous))
    )


@dataclass
class EnotsWolleyGenerator:
    """Transparent least-unused scan with fully pickleable mutable state."""

    terms: list[int] = field(default_factory=lambda: [1, 2])
    used: set[int] = field(default_factory=lambda: {1, 2})
    smallest_unused: int = 3

    def extend_to(self, count: int) -> None:
        if count < 0:
            raise ValueError("count must be nonnegative")
        if count <= len(self.terms):
            return
        if len(self.terms) < 2:
            raise RuntimeError("EnotsWolleyGenerator state is missing initial terms")

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


ENOTS_WOLLEY = SequenceDefinition[int](
    id="A336957",
    oeis="A336957",
    name="Enots--Wolley",
    aliases=("ew", "enots-wolley"),
    generator_factory=EnotsWolleyGenerator,
    generator_version=1,
    definition_version=1,
    offset=1,
    object_space=PositiveIntegers(),
    projections={"prime-exponents": prime_exponent_projection()},
    description=(
        "Lexicographically earliest unused positive-integer sequence starting 1, 2 "
        "where each later term shares a prime with its predecessor, is coprime to "
        "the term two places back, and introduces a prime absent from its predecessor."
    ),
)
