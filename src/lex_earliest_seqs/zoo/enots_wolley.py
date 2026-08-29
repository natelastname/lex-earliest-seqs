"""Enots--Wolley sequence definition (OEIS A336957)."""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cache
from heapq import heappop, heappush
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


def _next_coprime(lower_bound: int, forbidden_radical: int) -> int:
    """Return the least integer >= ``lower_bound`` coprime to ``forbidden_radical``."""

    candidate = max(1, lower_bound)
    while gcd(candidate, forbidden_radical) != 1:
        candidate += 1
    return candidate


@dataclass
class EnotsWolleyGenerator:
    """EW generator merging exact locally admissible candidate streams.

    For predecessor ``A`` and two-back term ``B``, let
    ``X = P(A) \\ P(B) = {p_1 < ... < p_k}``. Every integer that shares a prime
    with ``A`` and is coprime to ``B`` belongs to exactly one stream: stream ``i``
    consists of ``p_i * m`` with ``m`` coprime to ``rad(B) * p_1 * ... * p_{i-1}``.
    Thus the heap enumerates only values already satisfying the share/coprimality
    conditions. Used values and values introducing no new prime are skipped by
    advancing only their stream.

    Candidate support introduction is tested with a precomputed radical table.
    Used values remain a plain Python set. The radical table is derived and is
    therefore omitted from persisted pickles; it is rebuilt lazily only when a
    loaded generator is extended again. The persisted state layout is unchanged
    from generator version 3, so existing version-3 caches remain compatible.
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

    def _ensure_value(self, value: int) -> None:
        if value <= self.limit:
            self._ensure_radicals()
            return
        new_limit = self.limit
        while new_limit < value:
            new_limit *= 2
        self._resize(new_limit)

    def _next_candidate(self) -> int:
        previous = self.terms[-1]
        two_back = self.terms[-2]
        self._ensure_value(max(previous, two_back))
        assert self.radicals is not None

        previous_radical = self.radicals[previous]
        two_back_radical = self.radicals[two_back]
        shared_primes = tuple(sorted(prime_support(previous) - prime_support(two_back)))
        if not shared_primes:
            raise RuntimeError(
                "EW state has no predecessor prime disjoint from the two-back term"
            )

        # Each heap item is (candidate, stream_prime, multiplier,
        # forbidden_radical). Stream i is assigned exactly the candidates whose
        # least-indexed divisor from shared_primes is stream_prime.
        heap: list[tuple[int, int, int, int]] = []
        earlier_shared_product = 1
        for stream_prime in shared_primes:
            forbidden_radical = two_back_radical * earlier_shared_product
            lower_multiplier = (
                self.smallest_unused + stream_prime - 1
            ) // stream_prime
            multiplier = _next_coprime(lower_multiplier, forbidden_radical)
            heappush(
                heap,
                (
                    stream_prime * multiplier,
                    stream_prime,
                    multiplier,
                    forbidden_radical,
                ),
            )
            earlier_shared_product *= stream_prime

        while heap:
            candidate, stream_prime, multiplier, forbidden_radical = heappop(heap)
            self._ensure_value(candidate)
            assert self.radicals is not None

            # The stream construction already guarantees gcd(candidate, previous)
            # != 1 and gcd(candidate, two_back) == 1. The remaining greedy tests
            # are unusedness and introduction of a prime absent from previous.
            if (
                candidate not in self.used
                and previous_radical % self.radicals[candidate] != 0
            ):
                return candidate

            multiplier = _next_coprime(multiplier + 1, forbidden_radical)
            heappush(
                heap,
                (
                    stream_prime * multiplier,
                    stream_prime,
                    multiplier,
                    forbidden_radical,
                ),
            )

        raise RuntimeError("EW candidate heap unexpectedly exhausted")

    def extend_to(self, count: int) -> None:
        if count < 0:
            raise ValueError("count must be nonnegative")
        if count <= len(self.terms):
            return
        if len(self.terms) < 2:
            raise RuntimeError("EnotsWolleyGenerator state is missing initial terms")

        self._ensure_radicals()
        while len(self.terms) < count:
            candidate = self._next_candidate()
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
