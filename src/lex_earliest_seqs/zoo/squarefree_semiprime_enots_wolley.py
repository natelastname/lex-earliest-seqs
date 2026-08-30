"""Squarefree-semiprime Enots--Wolley sequence (OEIS placeholder A??????)."""

from __future__ import annotations

from dataclasses import dataclass, field

from ..core import SequenceDefinition
from ..object_space import PositiveIntegers
from ..projections import prime_exponent_projection, prime_factorization

# Direct greedy computation gives this finite transient.  A proved two-hub
# induction takes over afterward: for every k >= 10, the terms at indices
# 2k+5 and 2k+6 are 2*p_k and 3*p_k, with their order alternating with k.
_PROVED_PREFIX = (
    1,
    2,
    6,
    15,
    35,
    14,
    22,
    33,
    21,
    91,
    26,
    10,
    55,
    77,
    119,
    34,
    38,
    57,
    39,
    65,
    85,
    51,
    69,
    46,
)
_INITIAL_PRIMES = (2, 3, 5, 7, 11, 13, 17, 19, 23)


def is_squarefree_semiprime(value: int) -> bool:
    """Return whether ``value`` is a product of two distinct primes."""

    factors = prime_factorization(value)
    return len(factors) == 2 and all(exponent == 1 for _, exponent in factors)


@dataclass
class SquarefreeSemiprimeEnotsWolleyGenerator:
    """Generate the squarefree-semiprime restriction of Enots--Wolley.

    The initial terms are 1, 2.  Every later term is required to be a product
    of two distinct primes, while otherwise obeying the ordinary Enots--Wolley
    rule and lexicographic minimization.

    The greedy trajectory has a proved closed form after term 24.  If ``p_k``
    denotes the kth prime, then for every ``k >= 10``::

        {c_(2k+5), c_(2k+6)} = {2*p_k, 3*p_k},

    with ``(2*p_k, 3*p_k)`` for even k and the reverse order for odd k.  The
    generator therefore stores the finite greedy transient and uses the theorem
    for the infinite tail rather than repeatedly searching candidate integers.
    """

    terms: list[int] = field(default_factory=lambda: [1, 2])
    primes: list[int] = field(default_factory=lambda: list(_INITIAL_PRIMES))

    def _ensure_prime_count(self, count: int) -> None:
        if count <= len(self.primes):
            return

        candidate = self.primes[-1] + 2
        while len(self.primes) < count:
            is_prime = True
            for prime in self.primes:
                if prime * prime > candidate:
                    break
                if candidate % prime == 0:
                    is_prime = False
                    break
            if is_prime:
                self.primes.append(candidate)
            candidate += 2

    def _term_at_subscript(self, subscript: int) -> int:
        if subscript <= len(_PROVED_PREFIX):
            return _PROVED_PREFIX[subscript - 1]

        k = (subscript - 5) // 2
        self._ensure_prime_count(k)
        prime = self.primes[k - 1]
        first_subscript = 2 * k + 5
        first_factor = 2 if k % 2 == 0 else 3
        factor = first_factor if subscript == first_subscript else 5 - first_factor
        return factor * prime

    def extend_to(self, count: int) -> None:
        if count < 0:
            raise ValueError("count must be nonnegative")
        if count <= len(self.terms):
            return
        if len(self.terms) < 2:
            raise RuntimeError(
                "SquarefreeSemiprimeEnotsWolleyGenerator state is missing initial terms"
            )

        while len(self.terms) < count:
            subscript = len(self.terms) + 1
            self.terms.append(self._term_at_subscript(subscript))


SQUAREFREE_SEMIPRIME_ENOTS_WOLLEY = SequenceDefinition[int](
    id="squarefree-semiprime-enots-wolley",
    oeis="A??????",
    name="Squarefree-semiprime Enots--Wolley",
    aliases=(
        "semiprime-ew",
        "squarefree-semiprime-ew",
        "enots-wolley-squarefree-semiprime",
    ),
    generator_factory=SquarefreeSemiprimeEnotsWolleyGenerator,
    generator_version=1,
    definition_version=1,
    offset=1,
    object_space=PositiveIntegers(),
    projections={"prime-exponents": prime_exponent_projection()},
    description=(
        "Lexicographically earliest sequence starting 1, 2 and obeying the "
        "Enots--Wolley rule while requiring every later term to be a product "
        "of two distinct primes."
    ),
)
