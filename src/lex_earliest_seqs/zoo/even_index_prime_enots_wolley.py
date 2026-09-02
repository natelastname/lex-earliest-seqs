"""Enots--Wolley restricted to integers using an even-indexed prime."""

from __future__ import annotations

from dataclasses import dataclass, field
from math import isqrt

from ..core import SequenceDefinition
from ..object_space import PositiveIntegers
from ..projections import prime_exponent_projection
from .enots_wolley import is_candidate, prime_support
from .factor_restricted_enots_wolley import FactorRestrictedEnotsWolleyGenerator

# p_1 = 2, p_2 = 3, p_3 = 5, ... .  This table marks exactly the primes
# p_{2j}.  It is derived process-local state rather than generator state, so
# sequence caches do not pickle a potentially large sieve.
_even_index_prime_flags = bytearray(2)
_even_index_prime_limit = 1


def _ensure_even_index_prime_table(limit: int) -> None:
    """Ensure even-prime-index membership is tabulated through ``limit``."""

    global _even_index_prime_flags, _even_index_prime_limit

    if limit <= _even_index_prime_limit:
        return

    target = max(limit, 64 if _even_index_prime_limit < 64 else 2 * _even_index_prime_limit)
    sieve = bytearray(b"\x01") * (target + 1)
    sieve[0:2] = b"\x00\x00"
    for prime in range(2, isqrt(target) + 1):
        if not sieve[prime]:
            continue
        start = prime * prime
        sieve[start : target + 1 : prime] = b"\x00" * (
            (target - start) // prime + 1
        )

    flags = bytearray(target + 1)
    prime_index = 0
    for value in range(2, target + 1):
        if not sieve[value]:
            continue
        prime_index += 1
        if prime_index % 2 == 0:
            flags[value] = 1

    _even_index_prime_flags = flags
    _even_index_prime_limit = target


def is_even_index_prime(value: int) -> bool:
    """Return whether ``value`` is p_k for an even one-based prime index k."""

    if value < 2:
        return False
    _ensure_even_index_prime_table(value)
    return bool(_even_index_prime_flags[value])


@dataclass(frozen=True, slots=True)
class EvenIndexPrimePolicy:
    """Require a term to contain at least one even-indexed prime divisor."""

    def allows(self, value: int) -> bool:
        if value < 2:
            return False
        support = prime_support(value)
        if not support:
            return False
        _ensure_even_index_prime_table(max(support))
        return any(_even_index_prime_flags[prime] for prime in support)


@dataclass
class ReferenceEvenIndexPrimeEnotsWolleyGenerator:
    """Slow direct scanner used as a correctness oracle for X000006."""

    policy: EvenIndexPrimePolicy = field(default_factory=EvenIndexPrimePolicy)
    terms: list[int] = field(default_factory=lambda: [1, 3])
    used: set[int] = field(default_factory=lambda: {1, 3})

    def _next_candidate(self) -> int:
        previous = self.terms[-1]
        two_back = self.terms[-2]
        candidate = 2
        while True:
            if (
                candidate not in self.used
                and self.policy.allows(candidate)
                and is_candidate(candidate, previous, two_back)
            ):
                return candidate
            candidate += 1

    def extend_to(self, count: int) -> None:
        if count < 0:
            raise ValueError("count must be nonnegative")
        if count <= len(self.terms):
            return
        if len(self.terms) < 2:
            raise RuntimeError(
                "ReferenceEvenIndexPrimeEnotsWolleyGenerator state is missing "
                "initial terms"
            )

        while len(self.terms) < count:
            candidate = self._next_candidate()
            self.terms.append(candidate)
            self.used.add(candidate)


@dataclass
class EvenIndexPrimeEnotsWolleyGenerator(FactorRestrictedEnotsWolleyGenerator):
    """Optimized X000006 generator using persistent EW candidate streams.

    The inherited stream machinery only requires a global, history-independent
    ``policy.allows(value)`` predicate.  Here that predicate is membership in
    the integer universe having at least one divisor p_{2j}.  Odd-indexed primes
    remain legal as cofactors and may carry EW adjacency; they simply cannot make
    an integer admissible by themselves.
    """

    policy: EvenIndexPrimePolicy = field(default_factory=EvenIndexPrimePolicy)
    terms: list[int] = field(default_factory=lambda: [1, 3])
    used: set[int] = field(default_factory=lambda: {1, 3})


EVEN_INDEX_PRIME_ENOTS_WOLLEY = SequenceDefinition[int](
    id="X000006",
    oeis=None,
    name="Even-index-prime Enots--Wolley",
    aliases=(
        "even-index-prime-ew",
        "even-prime-index-ew",
        "alternating-prime-ew",
    ),
    generator_factory=EvenIndexPrimeEnotsWolleyGenerator,
    generator_version=1,
    definition_version=1,
    offset=1,
    object_space=PositiveIntegers(),
    projections={"prime-exponents": prime_exponent_projection()},
    description=(
        "Lexicographically earliest sequence starting 1, 3 and obeying the "
        "Enots--Wolley rule, with every later term required to be divisible by "
        "at least one even-indexed prime p_2, p_4, p_6, ... . Odd-indexed primes "
        "may still occur as cofactors."
    ),
)
