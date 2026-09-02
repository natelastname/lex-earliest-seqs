"""Enots--Wolley restricted to integers using every k-th prime."""

from __future__ import annotations

from array import array
from dataclasses import dataclass, field
from functools import partial
from math import isqrt

from ..core import SequenceDefinition
from ..object_space import PositiveIntegers
from ..projections import prime_exponent_projection
from .enots_wolley import is_candidate, prime_support
from .factor_restricted_enots_wolley import FactorRestrictedEnotsWolleyGenerator

# p_1 = 2, p_2 = 3, p_3 = 5, ... .  The table stores the one-based prime
# index at prime entries and 0 at composites.  It is derived process-local state
# rather than generator state, so sequence caches do not pickle a large sieve.
_prime_index_table = array("I", [0, 0])
_prime_index_limit = 1
_indexed_primes: list[int] = []


def _validate_k(k: int) -> None:
    if type(k) is not int:
        raise TypeError("k must be an integer")
    if k < 1:
        raise ValueError("k must be positive")


def _ensure_prime_index_table(limit: int) -> None:
    """Ensure one-based prime indices are tabulated through ``limit``."""

    global _prime_index_table, _prime_index_limit, _indexed_primes

    if limit <= _prime_index_limit:
        return

    target = max(limit, 64 if _prime_index_limit < 64 else 2 * _prime_index_limit)
    sieve = bytearray(b"\x01") * (target + 1)
    sieve[0:2] = b"\x00\x00"
    for prime in range(2, isqrt(target) + 1):
        if not sieve[prime]:
            continue
        start = prime * prime
        sieve[start : target + 1 : prime] = b"\x00" * (
            (target - start) // prime + 1
        )

    indices = array("I", [0]) * (target + 1)
    primes: list[int] = []
    prime_index = 0
    for value in range(2, target + 1):
        if not sieve[value]:
            continue
        prime_index += 1
        indices[value] = prime_index
        primes.append(value)

    _prime_index_table = indices
    _prime_index_limit = target
    _indexed_primes = primes


def prime_index(value: int) -> int | None:
    """Return the one-based index of prime ``value``, or ``None`` if nonprime."""

    if value < 2:
        return None
    _ensure_prime_index_table(value)
    index = int(_prime_index_table[value])
    return index or None


def nth_prime(index: int) -> int:
    """Return p_index for a positive one-based prime index."""

    _validate_k(index)
    target = max(64, _prime_index_limit)
    while len(_indexed_primes) < index:
        _ensure_prime_index_table(target)
        target *= 2
    return _indexed_primes[index - 1]


def is_every_kth_prime(value: int, k: int) -> bool:
    """Return whether ``value`` is p_j with j divisible by ``k``."""

    _validate_k(k)
    index = prime_index(value)
    return index is not None and index % k == 0


@dataclass(frozen=True, slots=True)
class EveryKthPrimePolicy:
    """Require a term to contain a prime p_j whose index j is a multiple of k."""

    k: int

    def __post_init__(self) -> None:
        _validate_k(self.k)

    def allows(self, value: int) -> bool:
        if value < 2:
            return False
        support = prime_support(value)
        if not support:
            return False
        _ensure_prime_index_table(max(support))
        return any(int(_prime_index_table[prime]) % self.k == 0 for prime in support)


@dataclass
class ReferenceEveryKthPrimeEnotsWolleyGenerator:
    """Slow direct scanner used as a correctness oracle for arbitrary k."""

    k: int = 2
    policy: EveryKthPrimePolicy = field(init=False)
    terms: list[int] = field(init=False)
    used: set[int] = field(init=False)

    def __post_init__(self) -> None:
        self.policy = EveryKthPrimePolicy(self.k)
        seed = nth_prime(self.k)
        self.terms = [1, seed]
        self.used = {1, seed}

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
                "ReferenceEveryKthPrimeEnotsWolleyGenerator state is missing "
                "initial terms"
            )

        while len(self.terms) < count:
            candidate = self._next_candidate()
            self.terms.append(candidate)
            self.used.add(candidate)


@dataclass
class EveryKthPrimeEnotsWolleyGenerator(FactorRestrictedEnotsWolleyGenerator):
    """Optimized every-k-th-prime EW generator using persistent candidate streams.

    ``k`` is arbitrary.  The allowed prime coordinates are
    ``p_k, p_{2k}, p_{3k}, ...``.  Every generated term must contain at least one
    such prime, while all other primes remain legal as cofactors and may carry
    EW adjacency.

    The inherited stream machinery only requires a global, history-independent
    ``policy.allows(value)`` predicate, so no k-specific candidate algorithm is
    needed.  Fresh instances seed themselves with ``1, p_k``.
    """

    policy: EveryKthPrimePolicy = field(init=False)
    terms: list[int] = field(init=False)
    used: set[int] = field(init=False)
    multiplier_successors: dict[int, dict[int, int]] = field(
        default_factory=dict,
        init=False,
        repr=False,
    )
    k: int = 2

    def __post_init__(self) -> None:
        self.policy = EveryKthPrimePolicy(self.k)
        seed = nth_prime(self.k)
        self.terms = [1, seed]
        self.used = {1, seed}


def make_every_kth_prime_enots_wolley_definition(
    *,
    id: str,
    k: int,
    name: str,
    aliases: tuple[str, ...] = (),
) -> SequenceDefinition[int]:
    """Build one registered member of the general every-k-th-prime EW family."""

    _validate_k(k)
    seed = nth_prime(k)
    return SequenceDefinition[int](
        id=id,
        oeis=None,
        name=name,
        aliases=aliases,
        generator_factory=partial(EveryKthPrimeEnotsWolleyGenerator, k=k),
        generator_version=1,
        definition_version=1,
        offset=1,
        object_space=PositiveIntegers(),
        projections={"prime-exponents": prime_exponent_projection()},
        description=(
            f"Lexicographically earliest sequence starting 1, {seed} and obeying "
            "the Enots--Wolley rule, with every later term required to be "
            f"divisible by at least one prime p_j whose one-based index j is a "
            f"multiple of {k}. Other primes remain legal as cofactors."
        ),
    )


EVERY_SECOND_PRIME_ENOTS_WOLLEY = make_every_kth_prime_enots_wolley_definition(
    id="X000006",
    k=2,
    name="Every-second-prime Enots--Wolley",
    aliases=(
        "every-second-prime-ew",
        "every-2nd-prime-ew",
        "even-index-prime-ew",
        "even-prime-index-ew",
        "alternating-prime-ew",
    ),
)

EVERY_THIRD_PRIME_ENOTS_WOLLEY = make_every_kth_prime_enots_wolley_definition(
    id="X000007",
    k=3,
    name="Every-third-prime Enots--Wolley",
    aliases=("every-third-prime-ew", "every-3rd-prime-ew"),
)

EVERY_FOURTH_PRIME_ENOTS_WOLLEY = make_every_kth_prime_enots_wolley_definition(
    id="X000008",
    k=4,
    name="Every-fourth-prime Enots--Wolley",
    aliases=("every-fourth-prime-ew", "every-4th-prime-ew"),
)

EVERY_KTH_PRIME_ENOTS_WOLLEY_DEFINITIONS = (
    EVERY_SECOND_PRIME_ENOTS_WOLLEY,
    EVERY_THIRD_PRIME_ENOTS_WOLLEY,
    EVERY_FOURTH_PRIME_ENOTS_WOLLEY,
)
