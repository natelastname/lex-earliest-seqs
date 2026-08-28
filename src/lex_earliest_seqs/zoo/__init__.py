"""Built-in sequence zoo."""

from __future__ import annotations

from ..core import SequenceRegistry
from .binary_enots_wolley import BINARY_ENOTS_WOLLEY, BinaryEnotsWolleyGenerator
from .enots_wolley import ENOTS_WOLLEY, EnotsWolleyGenerator

__all__ = [
    "BINARY_ENOTS_WOLLEY",
    "BinaryEnotsWolleyGenerator",
    "ENOTS_WOLLEY",
    "EnotsWolleyGenerator",
    "register_builtins",
]


def register_builtins(registry: SequenceRegistry) -> None:
    registry.register(ENOTS_WOLLEY)
    registry.register(BINARY_ENOTS_WOLLEY)
