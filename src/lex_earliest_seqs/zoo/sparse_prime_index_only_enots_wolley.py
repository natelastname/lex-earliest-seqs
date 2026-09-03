"""Enots--Wolley on sparse prime-index coordinate systems.

The allowed prime coordinates are selected by their one-based prime index.
Three registered families are provided:

- square indices: ``p_1, p_4, p_9, p_16, ...``;
- powers of two: ``p_1, p_2, p_4, p_8, ...``;
- self powers: ``p_1, p_4, p_27, p_256, ...``.

The production generator enumerates the retained multiplicative monoid by a
lazy heap merge.  It also carries prime-support bit masks for every multiplier
and selected term.  Candidate coprimality, novelty, and used-value retirement
therefore avoid factoring the often enormous sparse-coordinate terms.

The three index families share the candidate engine, but retain separate prime
frontiers.  In particular, a Rosser-style lower bound postpones exact ``p_n``
evaluation until an unopened coordinate can really beat the current monoid
head.  This is essential for the self-power family, whose next index jumps from
``7^7`` to ``8^8``.
"""

from __future__ import annotations

from bisect import bisect_left
from dataclasses import dataclass, field
from functools import partial
from heapq import heappop, heappush
from math import isqrt, log
from typing import Iterator, Literal

from ..core import SequenceDefinition
from ..object_space import PositiveIntegers
from ..projections import prime_exponent_projection
from .enots_wolley import is_candidate, prime_support
from .every_kth_prime_enots_wolley import nth_prime, prime_index

PrimeIndexFamily = Literal["square", "power_of_two", "self_power"]

_VALID_FAMILIES: frozenset[str] = frozenset(
    {"square", "power_of_two", "self_power"}
)


def _validate_family(family: str) -> None:
    if family not in _VALID_FAMILIES:
        allowed = ", ".join(sorted(_VALID_FAMILIES))
        raise ValueError(f"unknown prime-index family {family!r}; expected one of {allowed}")


def _iter_mask_positions(mask: int) -> Iterator[int]:
    """Yield set-bit positions of ``mask`` in increasing order."""

    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def _nth_prime_lower_bound(index: int) -> int:
    """Return a safe inexpensive lower bound for the one-based ``p_index``.

    For ``index >= 6`` we use a deliberately weakened form of the standard
    lower estimate ``p_n > n(log n + log log n - 1)``.  The ``-1.1`` margin is
    immaterial asymptotically but keeps this helper conservative.  A lower bound
    is sufficient: it lets the monoid merge defer an exact, potentially very
    expensive, nth-prime lookup while preserving sorted output.
    """

    if index < 1:
        raise ValueError("prime index must be positive")
    small = (2, 3, 5, 7, 11)
    if index <= len(small):
        return small[index - 1]
    n = float(index)
    estimate = int(n * (log(n) + log(log(n)) - 1.1))
    return max(index + 1, estimate)


def is_retained_prime_index(index: int, family: PrimeIndexFamily) -> bool:
    """Return whether one-based prime index ``index`` belongs to ``family``."""

    _validate_family(family)
    if type(index) is not int:
        raise TypeError("prime index must be an integer")
    if index < 1:
        return False

    if family == "square":
        root = isqrt(index)
        return root * root == index

    if family == "power_of_two":
        return index & (index - 1) == 0

    # self_power: 1^1, 2^2, 3^3, ...
    base = 1
    while True:
        value = base**base
        if value >= index:
            return value == index
        base += 1


def is_retained_prime(value: int, family: PrimeIndexFamily) -> bool:
    """Return whether prime ``value`` is retained by the selected index family."""

    index = prime_index(value)
    return index is not None and is_retained_prime_index(index, family)


@dataclass(frozen=True, slots=True)
class SparsePrimeIndexOnlyPolicy:
    """Allow exactly integers supported on one sparse prime-index family."""

    family: PrimeIndexFamily

    def __post_init__(self) -> None:
        _validate_family(self.family)

    def allows(self, value: int) -> bool:
        if value < 2:
            return False
        support = prime_support(value)
        return bool(support) and all(
            is_retained_prime(prime, self.family) for prime in support
        )


@dataclass
class ReferenceSparsePrimeIndexOnlyEnotsWolleyGenerator:
    """Very slow direct integer scanner used only as a tiny correctness oracle."""

    family: PrimeIndexFamily = "square"
    policy: SparsePrimeIndexOnlyPolicy = field(init=False)
    terms: list[int] = field(init=False)
    used: set[int] = field(init=False)

    def __post_init__(self) -> None:
        _validate_family(self.family)
        self.policy = SparsePrimeIndexOnlyPolicy(self.family)
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
        while len(self.terms) < count:
            candidate = self._next_candidate()
            self.terms.append(candidate)
            self.used.add(candidate)


@dataclass
class SparsePrimeIndexOnlyEnotsWolleyGenerator:
    """History-aware candidate merge for one sparse retained-prime monoid.

    ``multiplier_values`` is the increasing multiplicative monoid generated by
    the retained primes.  Its support masks use retained-prime *positions*, not
    ordinary prime indices: bit zero represents the first retained prime 2, bit
    one the next retained prime, and so on.

    The monoid heap merges the streams ``p * multiplier_values``.  Equal values
    from different prime representations are consumed together, so the shared
    table remains strictly increasing.  Candidate streams reuse that table and
    a per-retained-prime successor-with-delete structure for globally used
    products.
    """

    family: PrimeIndexFamily = "square"
    policy: SparsePrimeIndexOnlyPolicy = field(init=False)
    terms: list[int] = field(init=False)
    used: set[int] = field(init=False)
    term_support_masks: list[int] = field(init=False, repr=False)

    multiplier_values: list[int] = field(default_factory=lambda: [1], repr=False)
    multiplier_support_masks: list[int] = field(
        default_factory=lambda: [0],
        repr=False,
    )

    # Exact retained primes which have already been requested.  This cache may
    # contain one or more primes whose monoid streams are not opened yet.
    retained_primes: list[int] = field(default_factory=list, repr=False)

    # Heap entries are (product, retained_prime_position, multiplier_index).
    monoid_heap: list[tuple[int, int, int]] = field(default_factory=list, repr=False)
    next_allowed_position: int = 0

    # Per candidate-stream successor-with-delete tables over multiplier_values.
    # Keys are retained-prime positions, avoiding repeated large-prime hashing.
    multiplier_successors: dict[int, dict[int, int]] = field(
        default_factory=dict,
        repr=False,
    )

    def __post_init__(self) -> None:
        _validate_family(self.family)
        self.policy = SparsePrimeIndexOnlyPolicy(self.family)
        self.terms = [1, 2]
        self.used = {1, 2}
        self.term_support_masks = [0, 1]

    def _allowed_prime_index(self, position: int) -> int:
        """Return the ordinary prime index at retained position ``position``."""

        if position < 0:
            raise ValueError("retained-prime position must be nonnegative")
        base = position + 1
        if self.family == "square":
            return base * base
        if self.family == "power_of_two":
            return 1 << position
        return base**base

    def _prime_at_position(
        self,
        position: int,
        *,
        upper_bound: int | None = None,
    ) -> int | None:
        """Return a retained prime, optionally deferring an exact lookup.

        If ``upper_bound`` is supplied and a rigorous lower bound for the prime
        already exceeds it, ``None`` is returned without calling ``nth_prime``.
        Prime requests are otherwise filled in retained-position order.
        """

        if position < 0:
            raise ValueError("retained-prime position must be nonnegative")
        if position < len(self.retained_primes):
            return self.retained_primes[position]

        while len(self.retained_primes) <= position:
            next_position = len(self.retained_primes)
            index = self._allowed_prime_index(next_position)
            if (
                upper_bound is not None
                and next_position == position
                and _nth_prime_lower_bound(index) > upper_bound
            ):
                return None
            self.retained_primes.append(nth_prime(index))
        return self.retained_primes[position]

    def _peek_next_unopened_prime(self, upper_bound: int) -> int | None:
        return self._prime_at_position(
            self.next_allowed_position,
            upper_bound=upper_bound,
        )

    def _open_next_prime_stream(self) -> None:
        position = self.next_allowed_position
        prime = self._prime_at_position(position)
        assert prime is not None
        heappush(self.monoid_heap, (prime, position, 0))
        self.next_allowed_position += 1

    def _append_next_multiplier(self) -> None:
        """Append the next distinct retained-monoid multiplier and its mask."""

        while True:
            if not self.monoid_heap:
                self._open_next_prime_stream()
                continue

            open_min = self.monoid_heap[0][0]
            unopened = self._peek_next_unopened_prime(open_min)
            if unopened is not None and unopened <= open_min:
                self._open_next_prime_stream()
                continue

            value = open_min
            entries: list[tuple[int, int, int]] = []
            while self.monoid_heap and self.monoid_heap[0][0] == value:
                entries.append(heappop(self.monoid_heap))

            support_mask: int | None = None
            for _, prime_position, multiplier_index in entries:
                representation_mask = (
                    self.multiplier_support_masks[multiplier_index]
                    | (1 << prime_position)
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

            for _, prime_position, multiplier_index in entries:
                next_index = multiplier_index + 1
                if next_index >= len(self.multiplier_values):
                    raise RuntimeError("monoid stream advanced beyond shared table")
                prime = self.retained_primes[prime_position]
                heappush(
                    self.monoid_heap,
                    (
                        prime * self.multiplier_values[next_index],
                        prime_position,
                        next_index,
                    ),
                )
            return

    def _ensure_multiplier_index(self, index: int) -> None:
        if index < 0:
            raise ValueError("multiplier index must be nonnegative")
        while index >= len(self.multiplier_values):
            self._append_next_multiplier()

    def _multiplier_index_at_least(self, lower_bound: int) -> int:
        lower_bound = max(1, lower_bound)
        while self.multiplier_values[-1] < lower_bound:
            self._append_next_multiplier()
        return bisect_left(self.multiplier_values, lower_bound)

    def _find_multiplier_successor(self, stream_position: int, index: int) -> int:
        self._ensure_multiplier_index(index)
        parents = self.multiplier_successors.setdefault(stream_position, {})
        current = index
        path: list[int] = []
        while current in parents:
            path.append(current)
            current = parents[current]
            self._ensure_multiplier_index(current)
        for item in path:
            parents[item] = current
        return current

    def _delete_multiplier_index(self, stream_position: int, index: int) -> bool:
        surviving = self._find_multiplier_successor(stream_position, index)
        if surviving != index:
            return False
        parents = self.multiplier_successors.setdefault(stream_position, {})
        parents[index] = self._find_multiplier_successor(stream_position, index + 1)
        return True

    def _next_unused_multiplier_index(
        self,
        stream_position: int,
        index: int,
    ) -> int:
        parents = self.multiplier_successors.setdefault(stream_position, {})
        index = self._find_multiplier_successor(stream_position, index)
        stream_prime = self.retained_primes[stream_position]
        while True:
            multiplier = self.multiplier_values[index]
            if stream_prime * multiplier not in self.used:
                return index
            successor = self._find_multiplier_successor(stream_position, index + 1)
            parents[index] = successor
            index = successor

    def _next_stream_index(
        self,
        stream_position: int,
        lower_index: int,
        forbidden_mask: int,
    ) -> int:
        index = max(0, lower_index)
        while True:
            index = self._next_unused_multiplier_index(stream_position, index)
            if self.multiplier_support_masks[index] & forbidden_mask == 0:
                return index
            index += 1

    def _smallest_novel_prime(self, previous_mask: int, two_back_mask: int) -> int:
        """Return the least retained prime eligible to supply novelty."""

        blocked = previous_mask | two_back_mask
        position = 0
        while blocked & (1 << position):
            position += 1
        prime = self._prime_at_position(position)
        assert prime is not None
        return prime

    def _retire_used_value(self, value: int, support_mask: int) -> None:
        """Retire every already-materialized stream representation of ``value``."""

        table_max = self.multiplier_values[-1]
        for stream_position in _iter_mask_positions(support_mask):
            stream_prime = self.retained_primes[stream_position]
            multiplier = value // stream_prime
            if multiplier > table_max:
                continue
            index = bisect_left(self.multiplier_values, multiplier)
            if (
                index < len(self.multiplier_values)
                and self.multiplier_values[index] == multiplier
            ):
                self._delete_multiplier_index(stream_position, index)

    def _next_candidate(self) -> tuple[int, int]:
        previous_mask = self.term_support_masks[-1]
        two_back_mask = self.term_support_masks[-2]
        shared_mask = previous_mask & ~two_back_mask
        if shared_mask == 0:
            raise RuntimeError(
                "sparse prime-coordinate EW state has no predecessor prime "
                "disjoint from the two-back term"
            )

        # Every admissible candidate must introduce a retained prime absent from
        # both predecessor and two-back support.  No multiplier below the least
        # such prime can work, so all candidate streams can start at one shared
        # lower index instead of rescanning old-prime powers from 1.
        novel_prime = self._smallest_novel_prime(previous_mask, two_back_mask)
        lower_index = self._multiplier_index_at_least(novel_prime)

        # Each heap item is
        # (candidate, stream_position, multiplier_index, forbidden_mask).
        heap: list[tuple[int, int, int, int]] = []
        earlier_shared_mask = 0
        for stream_position in _iter_mask_positions(shared_mask):
            forbidden_mask = two_back_mask | earlier_shared_mask
            index = self._next_stream_index(
                stream_position,
                lower_index,
                forbidden_mask,
            )
            multiplier = self.multiplier_values[index]
            stream_prime = self.retained_primes[stream_position]
            heappush(
                heap,
                (
                    stream_prime * multiplier,
                    stream_position,
                    index,
                    forbidden_mask,
                ),
            )
            earlier_shared_mask |= 1 << stream_position

        while heap:
            candidate, stream_position, index, forbidden_mask = heappop(heap)
            multiplier_mask = self.multiplier_support_masks[index]
            candidate_mask = multiplier_mask | (1 << stream_position)

            # Stream construction guarantees allowed support, sharing with the
            # predecessor, two-back coprimality, unique stream ownership, and an
            # unused product.  Only predecessor-external novelty remains.
            if candidate_mask & ~previous_mask:
                return candidate, candidate_mask

            index = self._next_stream_index(
                stream_position,
                index + 1,
                forbidden_mask,
            )
            multiplier = self.multiplier_values[index]
            stream_prime = self.retained_primes[stream_position]
            heappush(
                heap,
                (
                    stream_prime * multiplier,
                    stream_position,
                    index,
                    forbidden_mask,
                ),
            )

        raise RuntimeError("sparse prime-coordinate candidate heap unexpectedly exhausted")

    def extend_to(self, count: int) -> None:
        if count < 0:
            raise ValueError("count must be nonnegative")
        if count <= len(self.terms):
            return

        while len(self.terms) < count:
            candidate, support_mask = self._next_candidate()
            self.terms.append(candidate)
            self.term_support_masks.append(support_mask)
            self.used.add(candidate)
            self._retire_used_value(candidate, support_mask)


def make_sparse_prime_index_only_enots_wolley_definition(
    *,
    id: str,
    family: PrimeIndexFamily,
    name: str,
    aliases: tuple[str, ...] = (),
    description: str,
) -> SequenceDefinition[int]:
    """Build one registered sparse-prime-index-only EW sequence."""

    _validate_family(family)
    return SequenceDefinition[int](
        id=id,
        oeis=None,
        name=name,
        aliases=aliases,
        generator_factory=partial(SparsePrimeIndexOnlyEnotsWolleyGenerator, family=family),
        generator_version=2,
        definition_version=1,
        offset=1,
        object_space=PositiveIntegers(),
        projections={"prime-exponents": prime_exponent_projection()},
        description=description,
    )


SQUARE_INDEX_PRIME_ONLY_ENOTS_WOLLEY = make_sparse_prime_index_only_enots_wolley_definition(
    id="X000012",
    family="square",
    name="Square-index-prime-only Enots--Wolley",
    aliases=(
        "square-index-prime-only-ew",
        "prime-coordinate-square-ew",
    ),
    description=(
        "Lexicographically earliest sequence starting 1, 2 and obeying the "
        "Enots--Wolley rule after retaining only prime coordinates "
        "p_1, p_4, p_9, p_16, ... ."
    ),
)

POWER_OF_TWO_INDEX_PRIME_ONLY_ENOTS_WOLLEY = (
    make_sparse_prime_index_only_enots_wolley_definition(
        id="X000013",
        family="power_of_two",
        name="Power-of-two-index-prime-only Enots--Wolley",
        aliases=(
            "power-of-two-index-prime-only-ew",
            "prime-coordinate-power-of-two-ew",
        ),
        description=(
            "Lexicographically earliest sequence starting 1, 2 and obeying the "
            "Enots--Wolley rule after retaining only prime coordinates "
            "p_1, p_2, p_4, p_8, p_16, ... ."
        ),
    )
)

SELF_POWER_INDEX_PRIME_ONLY_ENOTS_WOLLEY = make_sparse_prime_index_only_enots_wolley_definition(
    id="X000014",
    family="self_power",
    name="Self-power-index-prime-only Enots--Wolley",
    aliases=(
        "self-power-index-prime-only-ew",
        "prime-coordinate-self-power-ew",
    ),
    description=(
        "Lexicographically earliest sequence starting 1, 2 and obeying the "
        "Enots--Wolley rule after retaining only prime coordinates "
        "p_1, p_4, p_27, p_256, ... ."
    ),
)

SPARSE_PRIME_INDEX_ONLY_ENOTS_WOLLEY_DEFINITIONS = (
    SQUARE_INDEX_PRIME_ONLY_ENOTS_WOLLEY,
    POWER_OF_TWO_INDEX_PRIME_ONLY_ENOTS_WOLLEY,
    SELF_POWER_INDEX_PRIME_ONLY_ENOTS_WOLLEY,
)
