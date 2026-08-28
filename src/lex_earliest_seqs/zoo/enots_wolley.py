"""Enots--Wolley sequence definition (OEIS A336957)."""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cache
from math import gcd, isqrt

from ..core import SequenceDefinition
from ..object_space import PositiveIntegers
from ..projections import prime_exponent_projection

_INITIAL_LIMIT = 1_024


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
    """Return whether ``value`` satisfies the three EW adjacency conditions."""

    return (
        value >= 1
        and gcd(value, previous) != 1
        and gcd(value, two_back) == 1
        and bool(prime_support(value) - prime_support(previous))
    )


def _radical_table(limit: int) -> list[int]:
    """Return squarefree radicals for every integer through ``limit``."""

    if limit < 1:
        raise ValueError("limit must be positive")

    smallest_factor = list(range(limit + 1))
    smallest_factor[1] = 1
    for prime in range(2, isqrt(limit) + 1):
        if smallest_factor[prime] != prime:
            continue
        for value in range(prime * prime, limit + 1, prime):
            if smallest_factor[value] == value:
                smallest_factor[value] = prime

    radicals = [1] * (limit + 1)
    for value in range(2, limit + 1):
        prime = smallest_factor[value]
        quotient = value // prime
        radicals[value] = (
            radicals[quotient]
            if quotient % prime == 0
            else radicals[quotient] * prime
        )
    return radicals


@dataclass
class EnotsWolleyGenerator:
    """Fast least-unused EW scan with pickleable continuation state.

    Candidate support introduction is tested with a precomputed radical table.
    Used values remain a plain Python set. The radical table is derived and is
    therefore omitted from persisted pickles; it is rebuilt lazily only when a
    loaded generator is extended again.
    """

    terms: list[int] = field(default_factory=lambda: [1, 2])
    used: set[int] = field(default_factory=lambda: {1, 2})
    smallest_unused: int = 3
    limit: int = _INITIAL_LIMIT
    radicals: list[int] | None = field(default=None, repr=False)

    def __getstate__(self) -> dict[str, object]:
        """Persist continuation state while dropping the derived radical table."""

        state = self.__dict__.copy()
        state["radicals"] = None
        return state

    def _ensure_radicals(self) -> None:
        if self.radicals is None or len(self.radicals) != self.limit + 1:
            self.radicals = _radical_table(self.limit)

    def _resize(self, new_limit: int) -> None:
        if new_limit <= self.limit:
            return
        self.limit = new_limit
        self.radicals = _radical_table(new_limit)

    def extend_to(self, count: int) -> None:
        if count < 0:
            raise ValueError("count must be nonnegative")
        if count <= len(self.terms):
            return
        if len(self.terms) < 2:
            raise RuntimeError("EnotsWolleyGenerator state is missing initial terms")

        self._ensure_radicals()
        assert self.radicals is not None

        while len(self.terms) < count:
            previous = self.terms[-1]
            two_back = self.terms[-2]
            if previous > self.limit:
                new_limit = self.limit
                while new_limit < previous:
                    new_limit *= 2
                self._resize(new_limit)
                assert self.radicals is not None
            previous_radical = self.radicals[previous]
            candidate = self.smallest_unused

            while True:
                if candidate > self.limit:
                    self._resize(self.limit * 2)
                    assert self.radicals is not None
                    previous_radical = self.radicals[previous]

                # rad(candidate) divides rad(previous) exactly when candidate
                # introduces no prime absent from the predecessor.
                if (
                    candidate not in self.used
                    and gcd(candidate, previous) != 1
                    and gcd(candidate, two_back) == 1
                    and previous_radical % self.radicals[candidate] != 0
                ):
                    break
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
    generator_version=3,
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
