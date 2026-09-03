"""Family-specialized candidate engines for sparse prime-coordinate EW.

The baseline sparse generator keeps one stream for each retained predecessor
prime and scans the global retained multiplicative monoid until a multiplier
both avoids the two-back support and contains novelty.

For very sparse prime alphabets that still does unnecessary work: novelty is
known to contain at least one retained prime q outside the predecessor. This
module makes that coordinate explicit and searches canonical pair streams

    p_retained * q_new * m.

The factor q_new guarantees novelty, so the residual multiplier m can start at
1. This avoids materializing the global monoid all the way up to q_new merely
to discover q_new inside a multiplier.

The three registered families use two monoid backends:

* square-index primes keep the heap merge, which scales better once many
  retained prime streams are open;
* power-of-two and self-power indices use a simple pointer merge. Their prime
  alphabets are so sparse that a linear scan across the few open streams is
  cheaper and avoids heap duplicate traffic.

The self-power family additionally benefits from the pair-frontier cutoff: an
unopened p_{i^i} is evaluated exactly only if its rigorous lower bound can beat
an already available candidate.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

from .sparse_prime_index_only_enots_wolley import (
    POWER_OF_TWO_INDEX_PRIME_ONLY_ENOTS_WOLLEY as BASE_POWER_OF_TWO_DEFINITION,
    SELF_POWER_INDEX_PRIME_ONLY_ENOTS_WOLLEY as BASE_SELF_POWER_DEFINITION,
    SQUARE_INDEX_PRIME_ONLY_ENOTS_WOLLEY as BASE_SQUARE_DEFINITION,
    PrimeIndexFamily,
    SparsePrimeIndexOnlyEnotsWolleyGenerator,
    _iter_mask_positions,
)

PairKey = tuple[int, int]
PairForbiddenKey = tuple[int, int, int]


@dataclass
class PairFrontierSparsePrimeIndexOnlyEnotsWolleyGenerator(
    SparsePrimeIndexOnlyEnotsWolleyGenerator
):
    """Sparse EW candidate engine with novelty factored into pair streams.

    One candidate is assigned canonically to

    * its least eligible predecessor-retained prime ``p``; and
    * its least newly introduced retained prime ``q``.

    The residual multiplier may contain any retained primes except the two-back
    support and earlier canonical coordinates. Pair-local successor tables
    permanently skip residual multipliers whose full product has already been
    selected.

    A second cache remembers the first residual multiplier not already ruled out
    by each dynamic forbidden-support mask. Repeated local support states can
    therefore resume at their old frontier instead of rescanning the same prefix.

    The inherited single-prime ``multiplier_successors`` table is deliberately
    unused. A pair stream discovers a globally used product lazily and retires
    only that pair's residual multiplier.
    """

    pair_multiplier_successors: dict[PairKey, dict[int, int]] = field(
        default_factory=dict,
        repr=False,
    )
    pair_forbidden_cursors: dict[PairForbiddenKey, int] = field(
        default_factory=dict,
        repr=False,
    )

    def _retire_used_value(self, value: int, support_mask: int) -> None:
        """Do nothing; pair streams retire used products lazily."""

        del value, support_mask

    def _find_pair_successor(self, pair: PairKey, index: int) -> int:
        self._ensure_multiplier_index(index)
        parents = self.pair_multiplier_successors.setdefault(pair, {})
        current = index
        path: list[int] = []
        while current in parents:
            path.append(current)
            current = parents[current]
            self._ensure_multiplier_index(current)
        for item in path:
            parents[item] = current
        return current

    def _delete_pair_multiplier_index(self, pair: PairKey, index: int) -> int:
        parents = self.pair_multiplier_successors.setdefault(pair, {})
        current = self._find_pair_successor(pair, index)
        successor = self._find_pair_successor(pair, current + 1)
        parents[current] = successor
        return successor

    def _next_pair_multiplier_index(
        self,
        *,
        pair: PairKey,
        base_product: int,
        lower_index: int,
        forbidden_mask: int,
        candidate_limit: int | None,
    ) -> int | None:
        """Return the first residual multiplier giving a valid unused product.

        ``pair_forbidden_cursors`` is monotone for one dynamic support state:
        indices below the cursor were either structurally forbidden or have been
        permanently deleted as globally used pair products. The returned unused
        head itself is retained as the cursor, because another pair may win the
        global arbitration and leave this product available next time.
        """

        cursor_key = (pair[0], pair[1], forbidden_mask)
        index = max(lower_index, self.pair_forbidden_cursors.get(cursor_key, 0))
        index = self._find_pair_successor(pair, index)

        while True:
            multiplier = self.multiplier_values[index]
            candidate = base_product * multiplier
            if candidate_limit is not None and candidate >= candidate_limit:
                self.pair_forbidden_cursors[cursor_key] = index
                return None

            if candidate in self.used:
                index = self._delete_pair_multiplier_index(pair, index)
                self.pair_forbidden_cursors[cursor_key] = index
                continue

            if self.multiplier_support_masks[index] & forbidden_mask == 0:
                self.pair_forbidden_cursors[cursor_key] = index
                return index

            index = self._find_pair_successor(pair, index + 1)
            self.pair_forbidden_cursors[cursor_key] = index

    @staticmethod
    def _least_zero_position(mask: int) -> int:
        bit = ~mask & (mask + 1)
        return bit.bit_length() - 1

    def _next_candidate(self) -> tuple[int, int]:
        previous_mask = self.term_support_masks[-1]
        two_back_mask = self.term_support_masks[-2]
        shared_mask = previous_mask & ~two_back_mask
        if shared_mask == 0:
            raise RuntimeError(
                "sparse prime-coordinate EW state has no predecessor prime "
                "disjoint from the two-back term"
            )

        blocked_novel_mask = previous_mask | two_back_mask
        first_novel_position = self._least_zero_position(blocked_novel_mask)
        first_novel_prime = self._prime_at_position(first_novel_position)
        assert first_novel_prime is not None

        best_value: int | None = None
        best_mask = 0
        earlier_shared_mask = 0

        for stream_position in _iter_mask_positions(shared_mask):
            stream_prime = self.retained_primes[stream_position]

            # Every candidate assigned to this retained stream contains at least
            # first_novel_prime. Once even that base cannot beat the current
            # winner, later (larger) retained streams cannot beat it either.
            if (
                best_value is not None
                and stream_prime * first_novel_prime >= best_value
            ):
                break

            earlier_novel_mask = 0
            novel_position = 0
            while True:
                novel_bit = 1 << novel_position
                if blocked_novel_mask & novel_bit:
                    novel_position += 1
                    continue

                prime_upper_bound = None
                if best_value is not None:
                    prime_upper_bound = (best_value - 1) // stream_prime
                    if prime_upper_bound < 2:
                        break

                novel_prime = self._prime_at_position(
                    novel_position,
                    upper_bound=prime_upper_bound,
                )
                if novel_prime is None:
                    break

                base_product = stream_prime * novel_prime
                if best_value is not None and base_product >= best_value:
                    break

                pair = (stream_position, novel_position)
                forbidden_mask = (
                    two_back_mask | earlier_shared_mask | earlier_novel_mask
                )
                multiplier_index = self._next_pair_multiplier_index(
                    pair=pair,
                    base_product=base_product,
                    lower_index=0,
                    forbidden_mask=forbidden_mask,
                    candidate_limit=best_value,
                )
                if multiplier_index is not None:
                    multiplier = self.multiplier_values[multiplier_index]
                    candidate = base_product * multiplier
                    candidate_mask = (
                        self.multiplier_support_masks[multiplier_index]
                        | (1 << stream_position)
                        | novel_bit
                    )
                    if best_value is None or candidate < best_value:
                        best_value = candidate
                        best_mask = candidate_mask

                # A later q-stream is canonical only when the residual multiplier
                # omits all earlier novel coordinates. This removes duplicate
                # pair representations without changing the candidate set.
                earlier_novel_mask |= novel_bit
                novel_position += 1

            earlier_shared_mask |= 1 << stream_position

        if best_value is None:
            raise RuntimeError("sparse prime-coordinate pair frontiers exhausted")
        return best_value, best_mask


@dataclass
class PointerMonoidPairFrontierSparsePrimeIndexOnlyEnotsWolleyGenerator(
    PairFrontierSparsePrimeIndexOnlyEnotsWolleyGenerator
):
    """Pair-frontier engine using a pointer merge for very sparse prime sets."""

    # One multiplier-table pointer per opened retained-prime stream.
    monoid_stream_indices: list[int] = field(default_factory=list, repr=False)

    def _open_next_prime_stream(self) -> None:
        position = self.next_allowed_position
        prime = self._prime_at_position(position)
        assert prime is not None
        if position != len(self.monoid_stream_indices):
            raise RuntimeError("pointer monoid stream positions lost synchronization")
        self.monoid_stream_indices.append(0)
        self.next_allowed_position += 1

    def _append_next_multiplier(self) -> None:
        """Append the next monoid value by scanning the few open prime streams."""

        while True:
            if not self.monoid_stream_indices:
                self._open_next_prime_stream()
                continue

            open_min: int | None = None
            for position, multiplier_index in enumerate(self.monoid_stream_indices):
                prime = self.retained_primes[position]
                candidate = prime * self.multiplier_values[multiplier_index]
                if open_min is None or candidate < open_min:
                    open_min = candidate
            assert open_min is not None

            unopened = self._peek_next_unopened_prime(open_min)
            if unopened is not None and unopened <= open_min:
                self._open_next_prime_stream()
                continue

            value = open_min
            producers: list[int] = []
            support_mask: int | None = None
            for position, multiplier_index in enumerate(self.monoid_stream_indices):
                prime = self.retained_primes[position]
                if prime * self.multiplier_values[multiplier_index] != value:
                    continue
                producers.append(position)
                representation_mask = (
                    self.multiplier_support_masks[multiplier_index]
                    | (1 << position)
                )
                if support_mask is None:
                    support_mask = representation_mask
                elif representation_mask != support_mask:
                    raise RuntimeError("inconsistent support masks for monoid value")

            assert support_mask is not None
            if value <= self.multiplier_values[-1]:
                raise RuntimeError("retained monoid failed to increase strictly")
            self.multiplier_values.append(value)
            self.multiplier_support_masks.append(support_mask)

            for position in producers:
                next_index = self.monoid_stream_indices[position] + 1
                if next_index >= len(self.multiplier_values):
                    raise RuntimeError("pointer monoid advanced beyond shared table")
                self.monoid_stream_indices[position] = next_index
            return


@dataclass
class SquareIndexPrimeOnlyEnotsWolleyGenerator(
    PairFrontierSparsePrimeIndexOnlyEnotsWolleyGenerator
):
    """Square-index generator: pair frontiers plus scalable heap monoid merge."""

    family: PrimeIndexFamily = field(default="square", init=False)


@dataclass
class PowerOfTwoIndexPrimeOnlyEnotsWolleyGenerator(
    PointerMonoidPairFrontierSparsePrimeIndexOnlyEnotsWolleyGenerator
):
    """Power-of-two-index generator: pair frontiers plus small pointer merge."""

    family: PrimeIndexFamily = field(default="power_of_two", init=False)


@dataclass
class SelfPowerIndexPrimeOnlyEnotsWolleyGenerator(
    PointerMonoidPairFrontierSparsePrimeIndexOnlyEnotsWolleyGenerator
):
    """Self-power-index generator with pair-frontier prime-opening cutoffs."""

    family: PrimeIndexFamily = field(default="self_power", init=False)


SQUARE_INDEX_PRIME_ONLY_ENOTS_WOLLEY = replace(
    BASE_SQUARE_DEFINITION,
    generator_factory=SquareIndexPrimeOnlyEnotsWolleyGenerator,
    generator_version=3,
)

POWER_OF_TWO_INDEX_PRIME_ONLY_ENOTS_WOLLEY = replace(
    BASE_POWER_OF_TWO_DEFINITION,
    generator_factory=PowerOfTwoIndexPrimeOnlyEnotsWolleyGenerator,
    generator_version=3,
)

SELF_POWER_INDEX_PRIME_ONLY_ENOTS_WOLLEY = replace(
    BASE_SELF_POWER_DEFINITION,
    generator_factory=SelfPowerIndexPrimeOnlyEnotsWolleyGenerator,
    generator_version=3,
)

SPARSE_PRIME_INDEX_ONLY_ENOTS_WOLLEY_DEFINITIONS = (
    SQUARE_INDEX_PRIME_ONLY_ENOTS_WOLLEY,
    POWER_OF_TWO_INDEX_PRIME_ONLY_ENOTS_WOLLEY,
    SELF_POWER_INDEX_PRIME_ONLY_ENOTS_WOLLEY,
)
