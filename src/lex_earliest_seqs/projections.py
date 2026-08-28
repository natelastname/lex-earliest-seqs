"""Reusable incidence projections for common integer representations."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from functools import cache
from typing import Hashable, TypeVar

from .incidence import IncidenceProjection

FeatureT = TypeVar("FeatureT", bound=Hashable)


@cache
def prime_factorization(value: int) -> tuple[tuple[int, int], ...]:
    """Return a positive integer as sorted ``(prime, exponent)`` pairs."""

    if value < 1:
        raise ValueError(
            f"prime factorization requires a positive integer; got {value}"
        )

    factors: list[tuple[int, int]] = []
    remaining = value
    divisor = 2
    while divisor * divisor <= remaining:
        exponent = 0
        while remaining % divisor == 0:
            exponent += 1
            remaining //= divisor
        if exponent:
            factors.append((divisor, exponent))
        divisor = 3 if divisor == 2 else divisor + 2
    if remaining > 1:
        factors.append((remaining, 1))
    return tuple(factors)


def prime_coordinates(value: int) -> Mapping[int, int]:
    return dict(prime_factorization(value))


def primes_through(limit: int) -> tuple[int, ...]:
    if limit < 2:
        return ()
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for prime in range(2, int(limit**0.5) + 1):
        if not sieve[prime]:
            continue
        start = prime * prime
        sieve[start : limit + 1 : prime] = b"\x00" * (
            ((limit - start) // prime) + 1
        )
    return tuple(index for index, flag in enumerate(sieve) if flag)


def prime_exponent_projection() -> IncidenceProjection[int, int, int]:
    return IncidenceProjection(
        key="prime-exponents",
        title="Prime incidence chronology",
        coordinates=prime_coordinates,
        feature_sort_key=lambda prime: prime,
        through=lambda largest: primes_through(largest),
    )


def set_bit_positions(value: int) -> tuple[int, ...]:
    if value < 0:
        raise ValueError("binary digit projection requires a nonnegative integer")
    result: list[int] = []
    bit = 0
    remaining = value
    while remaining:
        if remaining & 1:
            result.append(bit)
        remaining >>= 1
        bit += 1
    return tuple(result)


def binary_digit_projection() -> IncidenceProjection[int, int, int]:
    return IncidenceProjection(
        key="binary-digits",
        title="Binary-digit incidence chronology",
        coordinates=lambda value: {bit: 1 for bit in set_bit_positions(value)},
        feature_sort_key=lambda bit: bit,
        feature_label=lambda bit: f"2^{bit}",
        through=lambda largest: range(largest + 1),
    )


def boolean_support_projection(
    *,
    key: str,
    title: str,
    support: Callable[[object], Iterable[FeatureT]],
    feature_sort_key: Callable[[FeatureT], object] | None = None,
) -> IncidenceProjection:
    sort_key = feature_sort_key or (lambda value: value)
    return IncidenceProjection(
        key=key,
        title=title,
        coordinates=lambda obj: {feature: 1 for feature in support(obj)},
        feature_sort_key=sort_key,
    )
