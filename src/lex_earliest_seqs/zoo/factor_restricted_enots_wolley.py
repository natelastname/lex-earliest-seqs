"""Reusable finite-Ω restrictions for Enots--Wolley-type sequences."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from functools import partial

from ..core import SequenceDefinition, SequenceGenerator
from ..object_space import PositiveIntegers
from ..projections import prime_exponent_projection, prime_factorization
from .enots_wolley import is_candidate


def big_omega(value: int) -> int:
    """Return Ω(value), counting prime factors with multiplicity."""

    if value < 1:
        raise ValueError("value must be positive")
    return sum(exponent for _, exponent in prime_factorization(value))


@dataclass(frozen=True, slots=True)
class EWFactorPolicy:
    """Finite multiplicative-degree restriction for an EW-type sequence.

    ``allowed_omega`` is a finite nonempty set of allowed values of the
    big-Omega function Ω(n), so multiplicity is counted.  ``squarefree`` adds
    the independent requirement that every prime exponent be at most one.

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
        if sum(exponent for _, exponent in factors) not in self.allowed_omega:
            return False
        return not self.squarefree or all(exponent == 1 for _, exponent in factors)


@dataclass
class FactorRestrictedEnotsWolleyGenerator:
    """Correct generic generator for a finite-Ω EW restriction.

    This is deliberately a simple reference/fallback implementation.  It scans
    admissible positive integers in numeric order at each step and applies the
    ordinary EW adjacency rule plus ``policy``.  Family members with exploitable
    structure may provide specialized generators while retaining the same
    ``EWFactorPolicy``.
    """

    policy: EWFactorPolicy
    terms: list[int] = field(default_factory=lambda: [1, 2])
    used: set[int] = field(default_factory=lambda: {1, 2})

    def _next_candidate(self) -> int:
        previous = self.terms[-1]
        two_back = self.terms[-2]
        candidate = 1
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
                "FactorRestrictedEnotsWolleyGenerator state is missing initial terms"
            )

        while len(self.terms) < count:
            candidate = self._next_candidate()
            self.terms.append(candidate)
            self.used.add(candidate)


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
) -> SequenceDefinition[int]:
    """Build metadata for one member of the finite-Ω EW family.

    When no specialized ``generator_factory`` is supplied, the generic direct
    generator is used.  Supplying an optimized generator changes only the
    implementation; ``policy`` remains the mathematical restriction defining
    the family member.
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
            f"Enots--Wolley rule with Ω(n) in {{{degrees}}}{squarefree}."
        )

    return SequenceDefinition[int](
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
    )
