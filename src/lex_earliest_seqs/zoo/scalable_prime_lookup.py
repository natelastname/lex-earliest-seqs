"""Exact nth-prime lookup for large sparse prime-coordinate queries.

Small prime indices are handled by the shared dense table in
:mod:`every_kth_prime_enots_wolley`, where the same sieve also supports inverse
prime-index queries. Large isolated indices are a different problem: sparse
coordinate families may ask directly for values such as ``p_(8**8)`` without
needing any intervening primes.

``primecountpy`` is therefore a required project dependency and the sole backend
here. Its ``nth_prime`` binding exposes primecount's direct nth-prime algorithm,
which combines prime counting near an approximation with a short final sieve.
There are deliberately no environment-sensitive executable/module fallbacks.
"""

from __future__ import annotations

from functools import cache

from primecountpy.primecount import nth_prime as primecount_nth_prime

# primecount's nth_prime() returns a signed 64-bit prime. This is exactly the
# number of primes below 2^63, and is enforced by primecount itself. Check it in
# Python first so callers get a stable, intelligible exception instead of a C++
# backend abort.
MAX_DIRECT_NTH_PRIME_INDEX = 216_289_611_853_439_384


def _validate_index(index: int) -> None:
    if type(index) is not int:
        raise TypeError("prime index must be an integer")
    if index < 1:
        raise ValueError("prime index must be positive")
    if index > MAX_DIRECT_NTH_PRIME_INDEX:
        raise OverflowError(
            "primecountpy nth_prime supports indices through "
            f"{MAX_DIRECT_NTH_PRIME_INDEX:,}; got {index:,}"
        )


@cache
def scalable_nth_prime(index: int) -> int:
    """Return the exact one-based ``index``-th prime via ``primecountpy``."""

    _validate_index(index)
    value = int(primecount_nth_prime(index))
    if value < 2:
        raise RuntimeError("primecountpy returned an invalid nth prime")
    return value


def _reset_for_tests() -> None:
    """Clear the process-local exact-query cache used by tests."""

    scalable_nth_prime.cache_clear()
