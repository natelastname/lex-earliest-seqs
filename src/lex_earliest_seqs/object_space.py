"""Common ordered ambient object spaces."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class PositiveIntegers:
    """Positive integers in their usual order: rank 0 is object 1."""

    key: str = "positive-integers"

    def at_rank(self, rank: int) -> int:
        if rank < 0:
            raise ValueError("rank must be nonnegative")
        return rank + 1

    def rank_of(self, value: int) -> int:
        if value < 1:
            raise ValueError("positive-integer objects must be at least 1")
        return value - 1


@dataclass(frozen=True, slots=True)
class NonNegativeIntegers:
    """Nonnegative integers in their usual order."""

    key: str = "nonnegative-integers"

    def at_rank(self, rank: int) -> int:
        if rank < 0:
            raise ValueError("rank must be nonnegative")
        return rank

    def rank_of(self, value: int) -> int:
        if value < 0:
            raise ValueError("nonnegative-integer objects must be nonnegative")
        return value
