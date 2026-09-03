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
from math import gcd, isqrt, prod

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


def _multiplier_introduces_new_prime(
    multiplier: int,
    previous_primes: frozenset[int],
) -> bool:
    """Test EW new-prime introduction without factoring the candidate.

    A candidate stream has the form ``r * multiplier`` where ``r`` is already
    in the predecessor support. Therefore the candidate introduces a new prime
    exactly when the multiplier has a prime factor outside that support. Strip
    all predecessor primes from the multiplier; a nontrivial remainder is
    equivalent to the required new prime.
    """

    remainder = multiplier
    for prime in previous_primes:
        while remainder % prime == 0:
            remainder //= prime
        if remainder == 1:
            return False
    return remainder > 1


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
    """EW streams with the full-return rule compiled into active streams.

    Two optimizations matter for this family.

    First, candidate new-prime testing never uses ordinary EW's dense radical
    table. For ``candidate = r*m`` the stream prime ``r`` is already in the
    predecessor support, so stripping predecessor primes from ``m`` is an exact
    new-prime test. This avoids maintaining a dense table out to potentially
    very large candidate values.

    Second, when the predecessor is target-free, an ordinary stream is split
    into two disjoint exact substreams: target-free candidates and full ``p*q``
    candidates. One-sided target candidates are therefore never put on the
    candidate heap. If the local lag-two/partition exclusions already contain
    either target prime, a full candidate is impossible and only the target-free
    substream is generated.

    Existing generator caches remain continuation-compatible: inherited
    successor maps still encode only permanently used products. Full substreams
    use their composite fixed factor ``r*p*q`` as an additional successor-map
    key, which cannot collide with the prime keys used by ordinary streams.
    """

    p: int = 2
    q: int = 3

    def __post_init__(self) -> None:
        self.p, self.q = _validate_pair(self.p, self.q)

    def extend_to(self, count: int) -> None:
        """Extend without allocating the ordinary EW radical table."""

        if count < 0:
            raise ValueError("count must be nonnegative")
        if count <= len(self.terms):
            return
        if len(self.terms) < 2:
            raise RuntimeError("FullReturnEnotsWolleyGenerator state is missing initial terms")

        while len(self.terms) < count:
            candidate = self._next_candidate()
            self.terms.append(candidate)
            self.used.add(candidate)
            while self.smallest_unused in self.used:
                self.smallest_unused += 1

    def _next_candidate(self) -> int:
        previous = self.terms[-1]
        two_back = self.terms[-2]

        previous_primes = prime_support(previous)
        two_back_primes = prime_support(two_back)
        two_back_radical = prod(two_back_primes)
        shared_primes = tuple(sorted(previous_primes - two_back_primes))
        if not shared_primes:
            raise RuntimeError(
                "full-return EW state has no predecessor prime disjoint from the two-back term"
            )

        force_full_return = full_return_restriction_active(previous, self.p, self.q)
        target_product = self.p * self.q

        # Heap item:
        #   (candidate, fixed_factor, multiplier, forbidden_radical,
        #    automatic_new_prime)
        #
        # Ordinary/target-free streams have fixed_factor equal to the shared
        # stream prime. A full-return substream has fixed_factor=r*p*q and its
        # multiplier is the remaining cofactor. Since p and q are absent from
        # a target-free predecessor, every such full candidate automatically
        # introduces new primes.
        heap: list[tuple[int, int, int, int, bool]] = []

        def push_stream(
            fixed_factor: int,
            lower_bound: int,
            forbidden_radical: int,
            automatic_new_prime: bool,
        ) -> None:
            multiplier = self._next_stream_multiplier(
                fixed_factor,
                lower_bound,
                forbidden_radical,
            )
            heappush(
                heap,
                (
                    fixed_factor * multiplier,
                    fixed_factor,
                    multiplier,
                    forbidden_radical,
                    automatic_new_prime,
                ),
            )

        earlier_shared_product = 1
        for stream_prime in shared_primes:
            forbidden_radical = two_back_radical * earlier_shared_product

            if not force_full_return:
                lower_multiplier = max(
                    1,
                    (self.smallest_unused + stream_prime - 1) // stream_prime,
                )
                push_stream(
                    stream_prime,
                    lower_multiplier,
                    forbidden_radical,
                    False,
                )
            else:
                # Target-free branch. Folding p*q into the coprimality
                # exclusion enumerates exactly the candidates containing neither
                # target prime. Duplicated factors are harmless to gcd().
                free_forbidden = forbidden_radical * target_product
                lower_multiplier = max(
                    1,
                    (self.smallest_unused + stream_prime - 1) // stream_prime,
                )
                push_stream(
                    stream_prime,
                    lower_multiplier,
                    free_forbidden,
                    False,
                )

                # Full-return branch. It exists only when lag-two and the
                # disjoint stream partition do not already forbid p or q.
                if gcd(forbidden_radical, target_product) == 1:
                    full_factor = stream_prime * target_product
                    lower_full_multiplier = max(
                        1,
                        (self.smallest_unused + full_factor - 1) // full_factor,
                    )
                    push_stream(
                        full_factor,
                        lower_full_multiplier,
                        forbidden_radical,
                        True,
                    )

            earlier_shared_product *= stream_prime

        while heap:
            (
                candidate,
                fixed_factor,
                multiplier,
                forbidden_radical,
                automatic_new_prime,
            ) = heappop(heap)

            introduces_new_prime = automatic_new_prime or _multiplier_introduces_new_prime(
                multiplier,
                previous_primes,
            )
            if introduces_new_prime:
                return candidate

            # _next_stream_multiplier already skips globally used products. A
            # new-prime failure is state-local, so advance only this heap stream
            # and do not install a permanent deletion for the rejected value.
            multiplier = self._next_stream_multiplier(
                fixed_factor,
                multiplier + 1,
                forbidden_radical,
            )
            heappush(
                heap,
                (
                    fixed_factor * multiplier,
                    fixed_factor,
                    multiplier,
                    forbidden_radical,
                    automatic_new_prime,
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
    id="X000015", p=2, q=3, aliases=("full-return-ew-2-3", "fr-ew-2-3")
)
FULL_RETURN_EW_2_5 = make_full_return_enots_wolley_definition(
    id="X000016", p=2, q=5, aliases=("full-return-ew-2-5", "fr-ew-2-5")
)
FULL_RETURN_EW_3_5 = make_full_return_enots_wolley_definition(
    id="X000017", p=3, q=5, aliases=("full-return-ew-3-5", "fr-ew-3-5")
)

FULL_RETURN_ENOTS_WOLLEY_DEFINITIONS = (
    FULL_RETURN_EW_2_3,
    FULL_RETURN_EW_2_5,
    FULL_RETURN_EW_3_5,
)
