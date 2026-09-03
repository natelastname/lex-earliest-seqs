"""Self-power-index-specific prime-coordinate optimization.

The retained prime indices for X000014 are fixed:

    1^1, 2^2, 3^3, 4^4, ...

Candidate arbitration sometimes has to inspect a future retained prime merely to
prove that it cannot beat the current best candidate.  For this family the
indices grow so quickly that such a losing probe can dominate the whole run:
without an optional primecount/primesieve backend, asking for p_(9^9) means a
portable segmented sieve out to the 387,420,489th prime.

The first ten self-power coordinates are tiny static mathematical data compared
with a persisted generator state.  Keeping them here makes those exact lookups
constant-time on every installation, while the generic scalable nth-prime backend
remains the fallback after the table is exhausted.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

from .scalable_prime_lookup import remember_nth_prime
from .sparse_prime_index_candidate_optimized import (
    SELF_POWER_INDEX_PRIME_ONLY_ENOTS_WOLLEY as BASE_DEFINITION,
    SelfPowerIndexPrimeOnlyEnotsWolleyGenerator as BaseSelfPowerIndexPrimeOnlyEnotsWolleyGenerator,
)
from .sparse_prime_index_only_enots_wolley import PrimeIndexFamily

# Exact p_(i^i) for i = 1,...,10.  These values were independently checked by
# verifying prime_pi(value) == i**i.  Position j corresponds to i = j + 1.
SELF_POWER_RETAINED_PRIME_PREFIX: tuple[int, ...] = (
    2,
    7,
    103,
    1_619,
    28_687,
    567_871,
    12_579_617,
    310_248_241,
    8_448_283_757,
    252_097_800_623,
)


@dataclass
class SelfPowerIndexPrimeOnlyEnotsWolleyGenerator(
    BaseSelfPowerIndexPrimeOnlyEnotsWolleyGenerator
):
    """Self-power EW with constant-time exact lookup for its early coordinates."""

    family: PrimeIndexFamily = field(default="self_power", init=False)

    def _prime_at_position(
        self,
        position: int,
        *,
        upper_bound: int | None = None,
    ) -> int | None:
        if position < 0:
            raise ValueError("retained-prime position must be nonnegative")
        if position < len(self.retained_primes):
            return self.retained_primes[position]

        if position >= len(SELF_POWER_RETAINED_PRIME_PREFIX):
            return super()._prime_at_position(position, upper_bound=upper_bound)

        # Because the table is exact, it gives a stronger cutoff than the generic
        # analytic lower bound: a future coordinate above the current admissible
        # upper bound can be rejected without materializing it at all.
        target_prime = SELF_POWER_RETAINED_PRIME_PREFIX[position]
        if upper_bound is not None and target_prime > upper_bound:
            return None

        while len(self.retained_primes) <= position:
            next_position = len(self.retained_primes)
            prime = SELF_POWER_RETAINED_PRIME_PREFIX[next_position]
            index = self._allowed_prime_index(next_position)
            remember_nth_prime(index, prime)
            self.retained_primes.append(prime)
        return self.retained_primes[position]


SELF_POWER_INDEX_PRIME_ONLY_ENOTS_WOLLEY = replace(
    BASE_DEFINITION,
    generator_factory=SelfPowerIndexPrimeOnlyEnotsWolleyGenerator,
    generator_version=4,
)
