"""Square-index-specific candidate optimization.

The pair-frontier candidate engine is a good fit for extremely sparse retained
prime alphabets, but square-index primes are already dense enough that opening a
separate ``(retained p, novel q)`` state is counterproductive.  At 10,000 terms
that representation created hundreds of thousands of dynamic pair states.

For the square family we keep the baseline single-retained-prime candidate
streams and the heap-merged retained monoid.  The only additional state is a
cursor for one exact local support configuration:

    (retained stream, predecessor support, forbidden support).

If the same local support state recurs, indices below the cursor have already
been proved structurally invalid (or permanently deleted as used products), so
candidate search resumes there instead of rescanning the same multiplier
prefix.  This preserves the compact candidate representation that empirically
wins for the square family.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

from .sparse_prime_index_only_enots_wolley import (
    SQUARE_INDEX_PRIME_ONLY_ENOTS_WOLLEY as BASE_SQUARE_DEFINITION,
    PrimeIndexFamily,
    SparsePrimeIndexOnlyEnotsWolleyGenerator,
)

StreamStateKey = tuple[int, int, int]


@dataclass
class SquareIndexPrimeOnlyEnotsWolleyGenerator(SparsePrimeIndexOnlyEnotsWolleyGenerator):
    """Square-index EW with cached local stream-scan frontiers."""

    family: PrimeIndexFamily = field(default="square", init=False)
    stream_state_cursors: dict[StreamStateKey, int] = field(
        default_factory=dict,
        repr=False,
    )

    def _next_valid_stream_index(
        self,
        stream_position: int,
        lower_index: int,
        forbidden_mask: int,
        previous_mask: int,
        candidate_limit: int | None,
    ) -> int | None:
        """Resume a repeated local support state at its proven scan frontier.

        The cache key contains every history-dependent predicate used below.
        Values below a saved cursor were either rejected because their support
        met ``forbidden_mask``, rejected because they supplied no novelty outside
        ``previous_mask``, or permanently deleted from this retained-prime stream
        after global selection.  None can become valid later for the same key.

        The current unused valid head itself remains the cursor.  If another
        stream wins the global least-choice arbitration, the head may still be
        needed when this support state recurs.
        """

        key = (stream_position, previous_mask, forbidden_mask)
        index = max(lower_index, self.stream_state_cursors.get(key, 0))
        stream_prime = self.retained_primes[stream_position]

        while True:
            index = self._next_unused_multiplier_index(stream_position, index)
            multiplier = self.multiplier_values[index]
            candidate = stream_prime * multiplier
            if candidate_limit is not None and candidate >= candidate_limit:
                self.stream_state_cursors[key] = index
                return None

            multiplier_mask = self.multiplier_support_masks[index]
            if (
                multiplier_mask & forbidden_mask == 0
                and multiplier_mask & ~previous_mask
            ):
                self.stream_state_cursors[key] = index
                return index

            index += 1
            self.stream_state_cursors[key] = index


SQUARE_INDEX_PRIME_ONLY_ENOTS_WOLLEY = replace(
    BASE_SQUARE_DEFINITION,
    generator_factory=SquareIndexPrimeOnlyEnotsWolleyGenerator,
    generator_version=3,
)
