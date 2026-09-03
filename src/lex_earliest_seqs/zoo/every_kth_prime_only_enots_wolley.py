"""Enots--Wolley on a thinned prime coordinate system.

For stride ``k``, retain the prime coordinates

    p_1, p_{1+k}, p_{1+2k}, ...

and forbid every other prime factor entirely.  Thus ``k=3`` retains
``2, 7, 17, 29, ...`` and no generated noninitial term is divisible by 3 or 5.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import partial

from ..core import SequenceDefinition
from ..object_space import PositiveIntegers
from ..projections import prime_exponent_projection
from .enots_wolley import is_candidate, prime_support
from .every_kth_prime_enots_wolley import prime_index
from .factor_restricted_enots_wolley import FactorRestrictedEnotsWolleyGenerator


def _validate_k(k: int) -> None:
    if type(k) is not int:
        raise TypeError("k must be an integer")
    if k < 1:
        raise ValueError("k must be positive")


def is_retained_prime(value: int, k: int) -> bool:
    """Return whether ``value`` is one of p_1, p_{1+k}, p_{1+2k}, ... ."""

    _validate_k(k)
    index = prime_index(value)
    return index is not None and (index - 1) % k == 0


@dataclass(frozen=True, slots=True)
class EveryKthPrimeOnlyPolicy:
    """Allow exactly integers whose prime support uses retained coordinates."""

    k: int

    def __post_init__(self) -> None:
        _validate_k(self.k)

    def allows(self, value: int) -> bool:
        if value < 2:
            return False
        support = prime_support(value)
        return bool(support) and all(is_retained_prime(prime, self.k) for prime in support)


@dataclass
class ReferenceEveryKthPrimeOnlyEnotsWolleyGenerator:
    """Slow direct scanner used as a correctness oracle for arbitrary k."""

    k: int = 2
    policy: EveryKthPrimeOnlyPolicy = field(init=False)
    terms: list[int] = field(init=False)
    used: set[int] = field(init=False)

    def __post_init__(self) -> None:
        self.policy = EveryKthPrimeOnlyPolicy(self.k)
        self.terms = [1, 2]
        self.used = {1, 2}

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
                "ReferenceEveryKthPrimeOnlyEnotsWolleyGenerator state is missing "
                "initial terms"
            )

        while len(self.terms) < count:
            candidate = self._next_candidate()
            self.terms.append(candidate)
            self.used.add(candidate)


@dataclass
class EveryKthPrimeOnlyEnotsWolleyGenerator(FactorRestrictedEnotsWolleyGenerator):
    """Optimized EW generator on retained prime coordinates only.

    ``k`` is arbitrary.  Fresh instances start ``1, 2``.  The inherited
    persistent candidate-stream machinery is valid because coordinate
    admissibility is a global, history-independent predicate: once a product
    contains a forbidden prime coordinate, that product can never become legal
    in a later EW state.
    """

    policy: EveryKthPrimeOnlyPolicy = field(init=False)
    terms: list[int] = field(init=False)
    used: set[int] = field(init=False)
    multiplier_successors: dict[int, dict[int, int]] = field(
        default_factory=dict,
        init=False,
        repr=False,
    )
    k: int = 2

    def __post_init__(self) -> None:
        self.policy = EveryKthPrimeOnlyPolicy(self.k)
        self.terms = [1, 2]
        self.used = {1, 2}


def make_every_kth_prime_only_enots_wolley_definition(
    *,
    id: str,
    k: int,
    name: str,
    aliases: tuple[str, ...] = (),
) -> SequenceDefinition[int]:
    """Build one registered member of the every-kth-prime-only EW family."""

    _validate_k(k)
    return SequenceDefinition[int](
        id=id,
        oeis=None,
        name=name,
        aliases=aliases,
        generator_factory=partial(EveryKthPrimeOnlyEnotsWolleyGenerator, k=k),
        generator_version=1,
        definition_version=1,
        offset=1,
        object_space=PositiveIntegers(),
        projections={"prime-exponents": prime_exponent_projection()},
        description=(
            "Lexicographically earliest sequence starting 1, 2 and obeying the "
            "Enots--Wolley rule after deleting all prime coordinates except "
            f"p_1, p_{{1+{k}}}, p_{{1+2*{k}}}, ... . Every prime factor of every "
            "later term must lie in that retained prime subsequence."
        ),
    )


EVERY_SECOND_PRIME_ONLY_ENOTS_WOLLEY = make_every_kth_prime_only_enots_wolley_definition(
    id="X000009",
    k=2,
    name="Every-second-prime-only Enots--Wolley",
    aliases=(
        "every-second-prime-only-ew",
        "prime-coordinate-ew-2",
        "prime-stride-ew-2",
    ),
)

EVERY_THIRD_PRIME_ONLY_ENOTS_WOLLEY = make_every_kth_prime_only_enots_wolley_definition(
    id="X000010",
    k=3,
    name="Every-third-prime-only Enots--Wolley",
    aliases=(
        "every-third-prime-only-ew",
        "prime-coordinate-ew-3",
        "prime-stride-ew-3",
    ),
)

EVERY_FOURTH_PRIME_ONLY_ENOTS_WOLLEY = make_every_kth_prime_only_enots_wolley_definition(
    id="X000011",
    k=4,
    name="Every-fourth-prime-only Enots--Wolley",
    aliases=(
        "every-fourth-prime-only-ew",
        "prime-coordinate-ew-4",
        "prime-stride-ew-4",
    ),
)

EVERY_KTH_PRIME_ONLY_ENOTS_WOLLEY_DEFINITIONS = (
    EVERY_SECOND_PRIME_ONLY_ENOTS_WOLLEY,
    EVERY_THIRD_PRIME_ONLY_ENOTS_WOLLEY,
    EVERY_FOURTH_PRIME_ONLY_ENOTS_WOLLEY,
)
