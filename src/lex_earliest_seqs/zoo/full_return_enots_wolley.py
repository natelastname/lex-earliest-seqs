"""Enots--Wolley with forced full returns for one target prime pair.

For a fixed pair of distinct primes ``p, q``, ordinary EW is unchanged except
when the immediately previous term is free of ``p`` and ``q``. In that state a
candidate involving either target prime is eligible only if it contains both:

    previous free of p,q
        => candidate is target-free or divisible by p*q.

Thus every first target return after a free term is full. One-sided target terms
remain possible whenever the predecessor already contains ``p`` or ``q``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import partial
from heapq import heappop, heappush
from math import isqrt

from ..core import SequenceDefinition
from ..object_space import PositiveIntegers
from ..projections import prime_exponent_projection
from .enots_wolley import EnotsWolleyGenerator, is_candidate, prime_support


def _is_prime(value: int) -> bool:
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    for divisor in range(3, isqrt(value) + 1, 2):
        if value % divisor == 0:
            return False
    return True


def _validate_pair(p: int, q: int) -> tuple[int, int]:
    if type(p) is not int or type(q) is not int:
        raise TypeError("p and q must be integers")
    if p == q:
        raise ValueError("p and q must be distinct")
    if not _is_prime(p) or not _is_prime(q):
        raise ValueError("p and q must both be prime")
    return (p, q) if p < q else (q, p)


def target_free(value: int, p: int, q: int) -> bool:
    return value % p != 0 and value % q != 0


def full_return_restriction_active(previous: int, p: int, q: int) -> bool:
    return target_free(previous, p, q)


def full_return_candidate_allowed(value: int, previous: int, p: int, q: int) -> bool:
    if not full_return_restriction_active(previous, p, q):
        return True
    return (value % p == 0) == (value % q == 0)


@dataclass
class ReferenceFullReturnEnotsWolleyGenerator:
    p: int = 2
    q: int = 3
    terms: list[int] = field(default_factory=lambda: [1, 2])
    used: set[int] = field(default_factory=lambda: {1, 2})
    smallest_unused: int = 3

    def __post_init__(self) -> None:
        self.p, self.q = _validate_pair(self.p, self.q)

    def _next_candidate(self) -> int:
        previous = self.terms[-1]
        two_back = self.terms[-2]
        candidate = self.smallest_unused
        while True:
            if (
                candidate not in self.used
                and is_candidate(candidate, previous, two_back)
                and full_return_candidate_allowed(candidate, previous, self.p, self.q)
            ):
                return candidate
            candidate += 1

    def extend_to(self, count: int) -> None:
        if count < 0:
            raise ValueError("count must be nonnegative")
        if count <= len(self.terms):
            return
        while len(self.terms) < count:
            candidate = self._next_candidate()
            self.terms.append(candidate)
            self.used.add(candidate)
            while self.smallest_unused in self.used:
                self.smallest_unused += 1


@dataclass
class FullReturnEnotsWolleyGenerator(EnotsWolleyGenerator):
    """Ordinary optimized EW streams plus the state-local full-return filter."""

    p: int = 2
    q: int = 3

    def __post_init__(self) -> None:
        self.p, self.q = _validate_pair(self.p, self.q)

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
                "full-return EW state has no predecessor prime disjoint from the two-back term"
            )

        force_full_return = full_return_restriction_active(previous, self.p, self.q)

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

            introduces_new_prime = (
                candidate not in self.used
                and previous_radical % self.radicals[candidate] != 0
            )
            return_allowed = True
            if force_full_return:
                return_allowed = (candidate % self.p == 0) == (candidate % self.q == 0)

            if introduces_new_prime and return_allowed:
                return candidate

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

        raise RuntimeError("full-return EW candidate heap unexpectedly exhausted")


def make_full_return_enots_wolley_definition(
    *,
    id: str,
    p: int,
    q: int,
    name: str | None = None,
    aliases: tuple[str, ...] = (),
) -> SequenceDefinition[int]:
    p, q = _validate_pair(p, q)
    return SequenceDefinition[int](
        id=id,
        oeis=None,
        name=name or f"Full-return Enots--Wolley ({p},{q})",
        aliases=aliases,
        generator_factory=partial(FullReturnEnotsWolleyGenerator, p=p, q=q),
        generator_version=1,
        definition_version=2,
        offset=1,
        object_space=PositiveIntegers(),
        projections={"prime-exponents": prime_exponent_projection()},
        description=(
            "Ordinary lexicographically earliest Enots--Wolley starting 1, 2, "
            f"with target pair ({p},{q}); whenever the immediately previous term "
            f"is divisible by neither {p} nor {q}, a candidate containing exactly "
            f"one of {p},{q} is ineligible, so the next target return must contain "
            "both target primes."
        ),
    )


FULL_RETURN_EW_2_3 = make_full_return_enots_wolley_definition(
    id="X000012", p=2, q=3, aliases=("full-return-ew-2-3", "fr-ew-2-3")
)
FULL_RETURN_EW_2_5 = make_full_return_enots_wolley_definition(
    id="X000013", p=2, q=5, aliases=("full-return-ew-2-5", "fr-ew-2-5")
)
FULL_RETURN_EW_3_5 = make_full_return_enots_wolley_definition(
    id="X000014", p=3, q=5, aliases=("full-return-ew-3-5", "fr-ew-3-5")
)

FULL_RETURN_ENOTS_WOLLEY_DEFINITIONS = (
    FULL_RETURN_EW_2_3,
    FULL_RETURN_EW_2_5,
    FULL_RETURN_EW_3_5,
)
