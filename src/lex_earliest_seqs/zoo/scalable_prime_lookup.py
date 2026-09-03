"""Scalable exact nth-prime lookup for sparse prime-coordinate sequences.

The ordinary prime-index table is ideal for small values because it supports
both ``p_n`` and inverse prime-index queries.  It is a poor fit for indices such
as ``8**8``: tabulating every integer through ``p_n`` consumes far more memory
than the sparse sequence itself.

This module provides a narrow exact ``nth_prime`` service:

1. use the optional :mod:`primesieve` C++ binding when it is installed;
2. otherwise advance a bounded-memory segmented odd sieve.

Both backends exploit the increasing query order of sparse prime-coordinate
families.  The C++ backend asks for the index gap after the nearest known lower
prime instead of restarting at zero.  The fallback keeps a process-global
forward cursor.  Out-of-order fallback queries remain correct; they use a
temporary fresh cursor rather than corrupting the forward state.
"""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from functools import cache
from importlib import import_module
from math import isqrt, log
from os import environ
from threading import Lock
from types import ModuleType

# One byte represents one odd integer, so this scans roughly two million
# integers per segment while using one MiB of temporary marking storage.
_SEGMENT_ODD_COUNT = 1 << 20

_PRIMESIEVE_UNCHECKED = object()
_primesieve_module: ModuleType | None | object = _PRIMESIEVE_UNCHECKED

# Exact backend queries are retained in sorted index order.  primesieve supports
# nth_prime(gap, start), so every later sparse coordinate can begin after the
# closest cached coordinate instead of recounting from 2.
_backend_lock = Lock()
_backend_indices: list[int] = [1]
_backend_primes: list[int] = [2]

_base_prime_limit = 1
_base_primes: list[int] = []

_segment_lock = Lock()
_segment_count = 1  # pi(2)
_segment_scanned_through = 2
_segment_requested_results: dict[int, int] = {1: 2}


def _validate_index(index: int) -> None:
    if type(index) is not int:
        raise TypeError("prime index must be an integer")
    if index < 1:
        raise ValueError("prime index must be positive")


def _load_primesieve() -> ModuleType | None:
    """Load the optional C++ backend once.

    ``LEX_EARLIEST_SEQS_DISABLE_PRIMESIEVE=1`` is useful for exercising the
    portable fallback in tests and benchmarks.
    """

    global _primesieve_module

    if environ.get("LEX_EARLIEST_SEQS_DISABLE_PRIMESIEVE") == "1":
        return None
    if _primesieve_module is _PRIMESIEVE_UNCHECKED:
        try:
            _primesieve_module = import_module("primesieve")
        except (ImportError, OSError):
            _primesieve_module = None
    assert _primesieve_module is None or isinstance(_primesieve_module, ModuleType)
    return _primesieve_module


def _nth_prime_upper_bound(index: int) -> int:
    """Return a safe elementary upper bound for ``p_index``."""

    small = (2, 3, 5, 7, 11)
    if index <= len(small):
        return small[index - 1]
    n = float(index)
    # Rosser's p_n < n(log n + log log n) bound holds for n >= 6.  The small
    # additive margin protects the integer conversion without affecting scale.
    return int(n * (log(n) + log(log(n)))) + 16


def _primes_through(limit: int) -> list[int]:
    """Return all ordinary primes through ``limit``, caching the largest sieve."""

    global _base_prime_limit, _base_primes

    if limit <= _base_prime_limit:
        return _base_primes[: bisect_right(_base_primes, limit)]

    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    for prime in range(2, isqrt(limit) + 1):
        if not sieve[prime]:
            continue
        start = prime * prime
        sieve[start : limit + 1 : prime] = b"\x00" * (
            (limit - start) // prime + 1
        )

    _base_primes = [value for value in range(2, limit + 1) if sieve[value]]
    _base_prime_limit = limit
    return _base_primes.copy()


def _advance_segmented(
    target_index: int,
    *,
    prime_count: int,
    scanned_through: int,
) -> tuple[int, int, int]:
    """Advance a segmented sieve cursor to exactly ``target_index``.

    Returns ``(prime, new_count, new_scanned_through)``.  The cursor deliberately
    stops at the requested prime instead of counting the rest of its final
    segment, making it valid to resume from that exact value.
    """

    if target_index < prime_count:
        raise ValueError("target index precedes segmented cursor")
    if target_index == prime_count:
        return scanned_through, prime_count, scanned_through

    upper = _nth_prime_upper_bound(target_index)
    low = max(3, scanned_through + 1)
    if low % 2 == 0:
        low += 1

    while True:
        base_primes = _primes_through(isqrt(upper) + 1)

        while low <= upper:
            high = min(upper, low + 2 * (_SEGMENT_ODD_COUNT - 1))
            if high % 2 == 0:
                high -= 1
            odd_count = (high - low) // 2 + 1
            flags = bytearray(b"\x01") * odd_count

            for prime in base_primes:
                if prime == 2:
                    continue
                start = max(prime * prime, ((low + prime - 1) // prime) * prime)
                if start % 2 == 0:
                    start += prime
                if start > high:
                    continue
                offset = (start - low) // 2
                flags[offset::prime] = b"\x00" * (
                    (odd_count - 1 - offset) // prime + 1
                )

            count_here = flags.count(1)
            if prime_count + count_here >= target_index:
                remaining = target_index - prime_count
                offset = -1
                for _ in range(remaining):
                    offset = flags.find(b"\x01", offset + 1)
                    if offset < 0:
                        raise RuntimeError("segmented prime count disagreed with flags")
                prime = low + 2 * offset
                return prime, target_index, prime

            prime_count += count_here
            scanned_through = high
            low = high + 2

        # The analytic upper bound should make this branch unreachable, but
        # retaining it makes correctness independent of a transcription error in
        # the bound or future changes to its cutoff.
        upper *= 2


def _segmented_nth_prime(index: int) -> int:
    global _segment_count, _segment_scanned_through

    with _segment_lock:
        cached = _segment_requested_results.get(index)
        if cached is not None:
            return cached

        if index < _segment_count:
            prime, _, _ = _advance_segmented(
                index,
                prime_count=1,
                scanned_through=2,
            )
            _segment_requested_results[index] = prime
            return prime

        prime, _segment_count, _segment_scanned_through = _advance_segmented(
            index,
            prime_count=_segment_count,
            scanned_through=_segment_scanned_through,
        )
        _segment_requested_results[index] = prime
        return prime


def _primesieve_nth_prime(backend: ModuleType, index: int) -> int:
    """Use the nearest known lower prime as primesieve's starting point."""

    with _backend_lock:
        insertion = bisect_left(_backend_indices, index)
        if insertion < len(_backend_indices) and _backend_indices[insertion] == index:
            return _backend_primes[insertion]

        if insertion == 0:
            value = int(backend.nth_prime(index))
        else:
            lower_index = _backend_indices[insertion - 1]
            lower_prime = _backend_primes[insertion - 1]
            value = int(backend.nth_prime(index - lower_index, lower_prime))

        if value < 2:
            raise RuntimeError("primesieve returned an invalid nth prime")
        _backend_indices.insert(insertion, index)
        _backend_primes.insert(insertion, value)
        return value


@cache
def scalable_nth_prime(index: int) -> int:
    """Return the one-based ``index``-th prime with bounded fallback memory."""

    _validate_index(index)
    backend = _load_primesieve()
    if backend is not None:
        return _primesieve_nth_prime(backend, index)
    return _segmented_nth_prime(index)


def _reset_for_tests() -> None:
    """Reset module caches; intended only for deterministic fallback tests."""

    global _primesieve_module
    global _backend_indices, _backend_primes
    global _base_prime_limit, _base_primes
    global _segment_count, _segment_scanned_through, _segment_requested_results

    with _backend_lock:
        _primesieve_module = _PRIMESIEVE_UNCHECKED
        _backend_indices = [1]
        _backend_primes = [2]
    with _segment_lock:
        _base_prime_limit = 1
        _base_primes = []
        _segment_count = 1
        _segment_scanned_through = 2
        _segment_requested_results = {1: 2}
        scalable_nth_prime.cache_clear()
