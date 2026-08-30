"""Reusable finite-omega restrictions for Enots--Wolley-type sequences."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from functools import partial
from heapq import heappop, heappush
from math import gcd, prod

from ..core import SequenceDefinition, SequenceGenerator
from ..object_space import PositiveIntegers
from ..projections import prime_exponent_projection, prime_factorization
from .enots_wolley import prime_support


def omega(value: int) -> int:
    """Return omega(value), the number of distinct prime factors."""

    if value < 1:
        raise ValueError("value must be positive")
    return len(prime_factorization(value))


def big_omega(value: int) -> int:
    """Return big-Omega(value), counting prime factors with multiplicity."""

    if value < 1:
        raise ValueError("value must be positive")
    return sum(exponent for _, exponent in prime_factorization(value))


@dataclass(frozen=True, slots=True)
class EWFactorPolicy:
    """Finite prime-support-size restriction for an EW-type sequence.

    ``allowed_omega`` is a finite nonempty set of allowed values of the
    little-omega function omega(n), so it counts distinct prime factors and
    ignores multiplicity. By default arbitrary positive exponents are allowed.
    ``squarefree=True`` independently requires every prime exponent to equal 1.

    For example, ``allowed_omega={2}, squarefree=False`` gives the biprimary
    family with exact supports {p, q} and objects p^a q^b for a,b >= 1, whereas
    setting ``squarefree=True`` leaves only pq on each pair support.

    The policy applies to generated terms after the sequence's initial seed;
    the standard EW seed 1, 2 need not itself satisfy the policy.
    """

    allowed_omega: frozenset[int]
    squarefree: bool = False

    def __post_init__(self) -> None:
        allowed = frozenset(self.allowed_omega)
        if not allowed:
            raise ValueError("allowed_omega must be nonempty")
        if any(type(degree) is not int for degree in allowed):
            raise TypeError("allowed_omega must contain only integers")
        if any(degree < 0 for degree in allowed):
            raise ValueError("allowed_omega values must be nonnegative")
        object.__setattr__(self, "allowed_omega", allowed)

    def allows(self, value: int) -> bool:
        """Return whether ``value`` satisfies this factor restriction."""

        if value < 1:
            return False
        factors = prime_factorization(value)
        if len(factors) not in self.allowed_omega:
            return False
        return not self.squarefree or all(exponent == 1 for _, exponent in factors)


@dataclass(frozen=True)
class FactorRestrictedEnotsWolleyDefinition(SequenceDefinition[int]):
    """Sequence definition carrying its finite-omega mathematical policy."""

    factor_policy: EWFactorPolicy | None = None

    def __post_init__(self) -> None:
        super().__post_init__()
        if self.factor_policy is None:
            raise ValueError("factor_policy must be provided")


@dataclass
class ReferenceFactorRestrictedEnotsWolleyGenerator:
    """Slow direct scanner used as a correctness oracle.

    This implementation deliberately knows nothing about candidate streams. At
    each step it scans positive integers from 1 upward and applies the complete
    mathematical definition. Production family members use
    ``FactorRestrictedEnotsWolleyGenerator`` instead, while tests compare the
    optimized generator against this reference implementation.
    """

    policy: EWFactorPolicy
    terms: list[int] = field(default_factory=lambda: [1, 2])
    used: set[int] = field(default_factory=lambda: {1, 2})

    def _is_candidate(self, value: int) -> bool:
        previous = self.terms[-1]
        two_back = self.terms[-2]
        value_support = prime_support(value)
        previous_support = prime_support(previous)
        return (
            value not in self.used
            and self.policy.allows(value)
            and bool(value_support & previous_support)
            and not bool(value_support & prime_support(two_back))
            and bool(value_support - previous_support)
        )

    def _next_candidate(self) -> int:
        candidate = 1
        while not self._is_candidate(candidate):
            candidate += 1
        return candidate

    def extend_to(self, count: int) -> None:
        if count < 0:
            raise ValueError("count must be nonnegative")
        if count <= len(self.terms):
            return
        if len(self.terms) < 2:
            raise RuntimeError(
                "ReferenceFactorRestrictedEnotsWolleyGenerator state is missing "
                "initial terms"
            )

        while len(self.terms) < count:
            candidate = self._next_candidate()
            self.terms.append(candidate)
            self.used.add(candidate)


def _next_coprime(lower_bound: int, forbidden_radical: int) -> int:
    candidate = max(1, lower_bound)
    while gcd(candidate, forbidden_radical) != 1:
        candidate += 1
    return candidate


@dataclass
class FactorRestrictedEnotsWolleyGenerator:
    """Optimized finite-omega EW generator using persistent candidate streams.

    For predecessor ``A`` and two-back term ``B``, let
    ``R = P(A) - P(B) = {p_1 < ... < p_k}``. Every integer sharing an active
    predecessor prime while remaining coprime to ``B`` is assigned to exactly
    one stream: stream ``i`` contains ``p_i*m`` where ``m`` is coprime to
    ``rad(B) * p_1 * ... * p_{i-1}``. The local streams are merged by a min-heap,
    so candidate integers are examined directly in increasing numerical order.

    The factor policy is global, not state-dependent. For each possible stream
    prime we therefore persist a successor-with-delete map over multipliers.
    Multipliers whose products are already used *or can never satisfy the factor
    policy* are permanently deleted and path-compressed. Only failures of the
    local "introduce a new predecessor-external prime" condition must be revisited
    in later states.

    Used-value retirement is eager across stream representations. Whenever a
    term ``x`` is selected, every representation ``x = p * (x/p)`` for
    ``p in P(x)`` is immediately deleted from the persistent successor map for
    stream ``p``. Since the current family has at most three distinct prime
    factors, this costs only O(omega(x)) <= 3 successor deletions per selected
    term and prevents the same used integer from being rediscovered later
    through another active-prime stream.

    This is the lazy version of the support-queue idea. Exact-support queues are
    not materialized wholesale: the integer stream reaches the current head of a
    structurally admissible support before any later element of that support can
    win. This avoids enumerating the very large collection of degree-3 supports
    whose radicals happen to lie below a coarse numerical bound.
    """

    policy: EWFactorPolicy
    terms: list[int] = field(default_factory=lambda: [1, 2])
    used: set[int] = field(default_factory=lambda: {1, 2})
    multiplier_successors: dict[int, dict[int, int]] = field(
        default_factory=dict,
        repr=False,
    )

    def _find_multiplier_successor(self, stream_prime: int, multiplier: int) -> int:
        parents = self.multiplier_successors.setdefault(stream_prime, {})
        current = max(1, multiplier)
        path: list[int] = []
        while current in parents:
            path.append(current)
            current = parents[current]
        for item in path:
            parents[item] = current
        return current

    def _delete_multiplier(self, stream_prime: int, multiplier: int) -> bool:
        """Permanently remove one multiplier from a stream successor set.

        Returns ``True`` only when this call performs a new deletion. Deletion is
        represented by linking ``multiplier`` to the first surviving successor,
        so future finds jump over it with path compression.
        """

        if multiplier < 1:
            raise ValueError("multiplier must be positive")

        surviving = self._find_multiplier_successor(stream_prime, multiplier)
        if surviving != multiplier:
            return False

        parents = self.multiplier_successors.setdefault(stream_prime, {})
        parents[multiplier] = self._find_multiplier_successor(
            stream_prime,
            multiplier + 1,
        )
        return True

    def _retire_used_value(self, value: int) -> int:
        """Delete every stream representation of a newly used value.

        If ``value = p*m`` and ``p`` is any prime divisor of ``value``, then the
        product can never be selected again from stream ``p`` in any future EW
        state. Retiring all such ``m`` eagerly prevents cross-stream
        rediscovery. The return value is the number of newly installed
        deletions, useful for tests and instrumentation.
        """

        if value < 1:
            raise ValueError("value must be positive")

        deletions = 0
        for stream_prime in prime_support(value):
            deletions += self._delete_multiplier(
                stream_prime,
                value // stream_prime,
            )
        return deletions

    def _next_persistently_eligible_multiplier(
        self,
        stream_prime: int,
        lower_bound: int,
    ) -> int:
        """Skip products that can never be selected in any future local state."""

        parents = self.multiplier_successors.setdefault(stream_prime, {})
        multiplier = self._find_multiplier_successor(stream_prime, lower_bound)
        while True:
            candidate = stream_prime * multiplier
            if candidate not in self.used and self.policy.allows(candidate):
                return multiplier
            successor = self._find_multiplier_successor(
                stream_prime,
                multiplier + 1,
            )
            parents[multiplier] = successor
            multiplier = successor

    def _next_stream_multiplier(
        self,
        stream_prime: int,
        lower_bound: int,
        forbidden_radical: int,
    ) -> int:
        """Return the least persistent candidate satisfying today's exclusions."""

        multiplier = max(1, lower_bound)
        while True:
            multiplier = self._next_persistently_eligible_multiplier(
                stream_prime,
                multiplier,
            )
            coprime_multiplier = _next_coprime(multiplier, forbidden_radical)
            if coprime_multiplier == multiplier:
                return multiplier
            multiplier = coprime_multiplier

    def _next_candidate(self) -> int:
        previous = self.terms[-1]
        two_back = self.terms[-2]
        previous_support = prime_support(previous)
        two_back_support = prime_support(two_back)
        shared_primes = tuple(sorted(previous_support - two_back_support))
        if not shared_primes:
            raise RuntimeError(
                "factor-restricted EW state has no predecessor prime disjoint "
                "from the two-back term"
            )

        # Stream i owns exactly the locally share/coprime integers whose first
        # divisor among shared_primes is p_i. The heap therefore merges disjoint
        # streams in exact numerical order.
        heap: list[tuple[int, int, int, int]] = []
        earlier_shared_product = 1
        two_back_radical = prod(two_back_support)
        for stream_prime in shared_primes:
            forbidden_radical = two_back_radical * earlier_shared_product
            multiplier = self._next_stream_multiplier(
                stream_prime,
                1,
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

            # Stream construction guarantees: candidate is unused, satisfies the
            # fixed factor policy, shares with the predecessor, and is coprime to
            # the two-back term. Only the local new-prime condition remains.
            if prime_support(candidate) - previous_support:
                return candidate

            multiplier = self._next_stream_multiplier(
                stream_prime,
                multiplier + 1,
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

        raise RuntimeError("factor-restricted EW candidate heap unexpectedly exhausted")

    def extend_to(self, count: int) -> None:
        if count < 0:
            raise ValueError("count must be nonnegative")
        if count <= len(self.terms):
            return
        if len(self.terms) < 2:
            raise RuntimeError(
                "FactorRestrictedEnotsWolleyGenerator state is missing initial terms"
            )

        while len(self.terms) < count:
            candidate = self._next_candidate()
            self.terms.append(candidate)
            self.used.add(candidate)
            self._retire_used_value(candidate)


def make_factor_restricted_enots_wolley_definition(
    *,
    id: str,
    name: str,
    policy: EWFactorPolicy,
    aliases: tuple[str, ...] = (),
    oeis: str | None = None,
    generator_factory: Callable[[], SequenceGenerator[int]] | None = None,
    generator_version: int = 1,
    definition_version: int = 1,
    description: str = "",
) -> FactorRestrictedEnotsWolleyDefinition:
    """Build metadata for one member of the finite-omega EW family.

    When no specialized ``generator_factory`` is supplied, the optimized
    persistent-stream generator is used. Tests and experiments can instantiate
    ``ReferenceFactorRestrictedEnotsWolleyGenerator`` explicitly as an oracle.
    Supplying another generator changes only the implementation; ``policy``
    remains the mathematical restriction defining the family member.
    """

    factory = generator_factory or partial(
        FactorRestrictedEnotsWolleyGenerator,
        policy=policy,
    )
    if not description:
        degrees = ", ".join(str(degree) for degree in sorted(policy.allowed_omega))
        squarefree = " and requiring squarefree terms" if policy.squarefree else ""
        description = (
            "Lexicographically earliest sequence starting 1, 2 and obeying the "
            f"Enots--Wolley rule with omega(n) in {{{degrees}}}{squarefree}."
        )

    return FactorRestrictedEnotsWolleyDefinition(
        id=id,
        oeis=oeis,
        name=name,
        aliases=aliases,
        generator_factory=factory,
        generator_version=generator_version,
        definition_version=definition_version,
        offset=1,
        object_space=PositiveIntegers(),
        projections={"prime-exponents": prime_exponent_projection()},
        description=description,
        factor_policy=policy,
    )
