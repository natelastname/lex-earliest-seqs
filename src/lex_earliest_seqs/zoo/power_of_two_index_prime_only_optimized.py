"""Power-of-two-index-specific candidate optimization.

For the retained prime coordinates

    p_1, p_2, p_4, p_8, p_16, ...

the explicit pair-frontier decomposition

    p_retained * q_new * m

is much faster than scanning the whole retained monoid for novelty. At 10,000
terms on the self-hosted research runner, the pair-frontier engine is about nine
times faster than the generic sparse candidate engine.

The retained monoid itself is still dense enough at that scale that the heap
merge narrowly beats the pointer merge (about 1.9 percent in the measured
10,000-term comparison). This module registers that measured winner rather than
sharing the self-power family's pointer backend.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

from .sparse_prime_index_candidate_optimized import (
    PairFrontierSparsePrimeIndexOnlyEnotsWolleyGenerator,
)
from .sparse_prime_index_only_enots_wolley import (
    POWER_OF_TWO_INDEX_PRIME_ONLY_ENOTS_WOLLEY as BASE_DEFINITION,
    PrimeIndexFamily,
)


@dataclass
class PowerOfTwoIndexPrimeOnlyEnotsWolleyGenerator(
    PairFrontierSparsePrimeIndexOnlyEnotsWolleyGenerator
):
    """Power-of-two-index EW using pair frontiers over the heap monoid merge."""

    family: PrimeIndexFamily = field(default="power_of_two", init=False)


POWER_OF_TWO_INDEX_PRIME_ONLY_ENOTS_WOLLEY = replace(
    BASE_DEFINITION,
    generator_factory=PowerOfTwoIndexPrimeOnlyEnotsWolleyGenerator,
    generator_version=3,
)
