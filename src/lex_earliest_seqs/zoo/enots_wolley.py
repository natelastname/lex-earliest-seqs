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
    """EW generator merging history-aware locally admissible candidate streams.

    For predecessor ``A`` and two-back term ``B``, let
    ``X = P(A) \\ P(B) = {p_1 < ... < p_k}``. Every integer that shares a prime
    with ``A`` and is coprime to ``B`` belongs to exactly one stream: stream ``i``
    consists of ``p_i * m`` with ``m`` coprime to ``rad(B) * p_1 * ... * p_{i-1}``.

    Each stream prime also owns a lazy successor-with-delete structure over its
    multipliers. Once ``p*m`` is discovered to be globally used, multiplier ``m``
    is permanently linked to its next possible successor. Path compression makes
    future occurrences of the same prime jump over historically exhausted
    multipliers rather than rediscovering them. Local coprimality failures and
    failures to introduce a new prime are never deleted because those conditions
    can change with the predecessor pair.

    Candidate support introduction is tested with a precomputed radical table.
    Used values remain a plain Python set. The radical table is derived and is
    therefore omitted from persisted pickles; the multiplier-successor maps are
    persisted. Older generator-version-3 caches have no successor maps, so they
    acquire them lazily on first extension and remain compatible.
    """

    terms: list[int] = field(default_factory=lambda: [1, 2])
    used: set[int] = field(default_factory=lambda: {1, 2})
    smallest_unused: int = 3
    limit: int = _INITIAL_LIMIT
    radicals: list[int] | None = field(default=None, repr=False)
    unused_multiplier_successors: dict[int, dict[int, int]] = field(
        default_factory=dict,
        repr=False,
    )

    def __getstate__(self) -> dict[str, object]:
        """Persist continuation/history state while dropping derived radicals."""

        state = self.__dict__.copy()
        state["radicals"] = None
        return state

    def _successor_maps(self) -> dict[int, dict[int, int]]:
        """Return successor maps, initializing them for older version-3 pickles."""

        maps = getattr(self, "unused_multiplier_successors", None)
        if maps is None:
            maps = {}
            self.unused_multiplier_successors = maps
        return maps

    def _find_multiplier_successor(self, stream_prime: int, multiplier: int) -> int:
        """Find the first multiplier not lazily deleted for ``stream_prime``."""

        parents = self._successor_maps().setdefault(stream_prime, {})
        current = max(1, multiplier)
        path: list[int] = []
        while current in parents:
            path.append(current)
            current = parents[current]
        for item in path:
            parents[item] = current
        return current

    def _next_unused_multiplier(self, stream_prime: int, lower_bound: int) -> int:
        """Return least ``m >= lower_bound`` not known to have ``p*m`` used.

        Used products are permanent, so every discovered used multiplier is
        safely deleted from this prime's successor set forever.
        """

        parents = self._successor_maps().setdefault(stream_prime, {})
        multiplier = self._find_multiplier_successor(stream_prime, lower_bound)
        while stream_prime * multiplier in self.used:
            successor = self._find_multiplier_successor(
                stream_prime,
                multiplier + 1,
            )
            parents[multiplier] = successor
            multiplier = successor
        return multiplier

    def _next_stream_multiplier(
        self,
        stream_prime: int,
        lower_bound: int,
        forbidden_radical: int,
    ) -> int:
        """Return least multiplier satisfying global-unused and local-coprime tests."""

        multiplier = max(1, lower_bound)
        while True:
            multiplier = self._next_unused_multiplier(stream_prime, multiplier)
            coprime_multiplier = _next_coprime(multiplier, forbidden_radical)
            if coprime_multiplier == multiplier:
                return multiplier
            multiplier = coprime_multiplier

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
            lower_multiplier = max(
                1,
                (self.smallest_unused + stream_prime - 1) // stream_prime,
            )
            multiplier = self._next_stream_multiplier(
                stream_prime,
                lower_multiplier,
                forbidden_radical,
            )
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

            # Stream construction guarantees the share/two-back-coprimality
            # conditions and skips products already known to be used. The only
            # mathematical condition left is introduction of a new prime.
            if (
                candidate not in self.used
                and previous_radical % self.radicals[candidate] != 0
            ):
                return candidate

            # A used candidate may be permanently deleted by the successor DSU.
            # A new-prime failure must only be skipped for this local state.
            next_lower_bound = multiplier if candidate in self.used else multiplier + 1
            multiplier = self._next_stream_multiplier(
                stream_prime,
                next_lower_bound,
                forbidden_radical,
            )
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
