"""Forced-squarefree Enots--Wolley sequence (OEIS A399457)."""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cache

from ..core import SequenceDefinition
from ..object_space import PositiveIntegers
from ..projections import prime_exponent_projection, prime_factorization
from .enots_wolley import is_candidate


@cache
def is_squarefree(value: int) -> bool:
    """Return whether ``value`` is squarefree, with 1 squarefree by convention."""

    if value < 1:
        return False
    return all(exponent == 1 for _, exponent in prime_factorization(value))


@dataclass
class ForcedSquarefreeEnotsWolleyGenerator:
    """Transparent EW scan permanently excluding nonsquarefree candidates."""

    terms: list[int] = field(default_factory=lambda: [1, 2])
    used: set[int] = field(default_factory=lambda: {1, 2})
    smallest_unused_squarefree: int = 3

    def extend_to(self, count: int) -> None:
        if count < 0:
            raise ValueError("count must be nonnegative")
        if count <= len(self.terms):
            return
        if len(self.terms) < 2:
            raise RuntimeError(
                "ForcedSquarefreeEnotsWolleyGenerator state is missing initial terms"
            )

        while len(self.terms) < count:
            candidate = self.smallest_unused_squarefree
            while (
                candidate in self.used
                or not is_squarefree(candidate)
                or not is_candidate(candidate, self.terms[-1], self.terms[-2])
            ):
                candidate += 1

            self.terms.append(candidate)
            self.used.add(candidate)

            while (
                self.smallest_unused_squarefree in self.used
                or not is_squarefree(self.smallest_unused_squarefree)
            ):
                self.smallest_unused_squarefree += 1


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
    generator_version=1,
    definition_version=1,
    offset=1,
    object_space=PositiveIntegers(),
    projections={"prime-exponents": prime_exponent_projection()},
    description=(
        "Lexicographically earliest sequence obeying the Enots--Wolley rule while "
        "permanently excluding nonsquarefree positive integers."
    ),
)
