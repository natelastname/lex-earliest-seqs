"""Forced-squarefree Enots--Wolley sequence (OEIS A399457)."""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cache
from heapq import heappop, heappush
from math import gcd, prod

from ..core import SequenceDefinition
from ..object_space import PositiveIntegers
from ..projections import prime_exponent_projection, prime_factorization
from .enots_wolley import _radical_table, prime_support

_INITIAL_LIMIT = 1_024


@cache
def is_squarefree(value: int) -> bool:
    """Return whether ``value`` is squarefree, with 1 squarefree by convention."""

    if value < 1:
        return False
    return all(exponent == 1 for _, exponent in prime_factorization(value))


def _nonempty_subset_products(primes: tuple[int, ...]) -> list[int]:
    """Return products of all nonempty subsets of distinct ``primes``."""

    products = [1]
    for prime in primes:
        products += [value * prime for value in products]
    return products[1:]


@dataclass
class ForcedSquarefreeEnotsWolleyGenerator:
    """Forced-squarefree EW generator with history-aware candidate streams.

    For predecessor A and two-back term B, every admissible squarefree candidate
    has a unique factorization ``x*y``. ``x`` is a nonempty squarefree product of
    primes in P(A) - P(B), while ``y > 1`` is squarefree and coprime to rad(A*B).
    Each possible ``x`` defines one increasing candidate stream.

    ``cofactor_frontiers[x]`` remembers the least base-eligible cofactor that has
    not already produced a used product ``x*y``. This frontier is independent of
    the current A/B pair, so repeated appearances of the same ``x`` never rescan
    an exhausted historical prefix. Current local coprimality restrictions are
    imposed only after jumping to that persistent frontier.

    The radical table accelerates squarefree/coprimality successors. It is
    derived state and is therefore omitted from persisted pickles; the history
    frontiers themselves are persisted because rebuilding them would discard the
    main late-term acceleration.
    """

    terms: list[int] = field(default_factory=lambda: [1, 2])
    used: set[int] = field(default_factory=lambda: {1, 2})
    cofactor_frontiers: dict[int, int] = field(default_factory=dict)
    limit: int = _INITIAL_LIMIT
    radicals: list[int] | None = field(default=None, repr=False)

    def __getstate__(self) -> dict[str, object]:
        state = self.__dict__.copy()
        state["radicals"] = None
        return state

    def _ensure_radicals(self) -> None:
        if self.radicals is None or len(self.radicals) != self.limit + 1:
            self.radicals = _radical_table(self.limit)

    def _ensure_value(self, value: int) -> None:
        if value <= self.limit:
            self._ensure_radicals()
            return
        new_limit = self.limit
        while new_limit < value:
            new_limit *= 2
        self.limit = new_limit
        self.radicals = _radical_table(new_limit)

    def _next_squarefree_coprime(self, lower_bound: int, forbidden_radical: int) -> int:
        """Return least y >= lower_bound squarefree and coprime to forbidden_radical."""

        candidate = max(2, lower_bound)
        while True:
            self._ensure_value(candidate)
            assert self.radicals is not None
            if (
                self.radicals[candidate] == candidate
                and gcd(candidate, forbidden_radical) == 1
            ):
                return candidate
            candidate += 1

    def _cofactor_frontier(self, shared_part: int) -> int:
        """Return the least globally unused base cofactor for ``shared_part``.

        The base stream for ``x`` consists of squarefree ``y > 1`` coprime to
        ``x``. Once ``x*y`` has appeared, that cofactor can never win again under
        any future local EW state, so the stored frontier may advance past it
        permanently. Local exclusions from the current predecessor/two-back pair
        are deliberately *not* folded into this persistent frontier.
        """

        outside_part = max(2, self.cofactor_frontiers.get(shared_part, 2))
        while True:
            outside_part = self._next_squarefree_coprime(
                outside_part,
                shared_part,
            )
            if shared_part * outside_part not in self.used:
                self.cofactor_frontiers[shared_part] = outside_part
                return outside_part
            outside_part += 1

    def _next_candidate(self) -> int:
        previous = self.terms[-1]
        two_back = self.terms[-2]
        previous_support = prime_support(previous)
        two_back_support = prime_support(two_back)
        shared_primes = tuple(sorted(previous_support - two_back_support))
        if not shared_primes:
            raise RuntimeError(
                "forced-squarefree EW state has no predecessor prime disjoint "
                "from the two-back term"
            )

        forbidden_radical = prod(previous_support | two_back_support)
        heap: list[tuple[int, int, int]] = []

        for shared_part in _nonempty_subset_products(shared_primes):
            # First jump over the entire globally exhausted prefix for this x.
            # Then impose today's stronger coprimality condition. A cofactor that
            # is blocked only by the current A/B pair remains available to future
            # states, so it does not advance the persistent frontier.
            outside_part = self._next_squarefree_coprime(
                self._cofactor_frontier(shared_part),
                forbidden_radical,
            )
            heappush(
                heap,
                (shared_part * outside_part, shared_part, outside_part),
            )

        while heap:
            candidate, shared_part, outside_part = heappop(heap)
            if candidate not in self.used:
                return candidate

            outside_part = self._next_squarefree_coprime(
                outside_part + 1,
                forbidden_radical,
            )
            heappush(
                heap,
                (shared_part * outside_part, shared_part, outside_part),
            )

        raise RuntimeError("forced-squarefree EW candidate heap unexpectedly exhausted")

    def extend_to(self, count: int) -> None:
        if count < 0:
            raise ValueError("count must be nonnegative")
        if count <= len(self.terms):
            return
        if len(self.terms) < 2:
            raise RuntimeError(
                "ForcedSquarefreeEnotsWolleyGenerator state is missing initial terms"
            )

        self._ensure_radicals()
        while len(self.terms) < count:
            candidate = self._next_candidate()
            self.terms.append(candidate)
            self.used.add(candidate)


FORCED_SQUAREFREE_ENOTS_WOLLEY = SequenceDefinition[int](
    id="A399457",
    oeis="A399457",
    name="Forced-squarefree Enots--Wolley",
    aliases=(
        "forced-squarefree-ew",
        "squarefree-ew",
        "enots-wolley-forced-squarefree",
    ),
    generator_factory=ForcedSquarefreeEnotsWolleyGenerator,
    generator_version=3,
    definition_version=1,
    offset=1,
    object_space=PositiveIntegers(),
    projections={"prime-exponents": prime_exponent_projection()},
    description=(
        "Lexicographically earliest sequence obeying the Enots--Wolley rule while "
        "permanently excluding nonsquarefree positive integers."
    ),
)
