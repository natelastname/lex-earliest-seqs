"""Built-in sequence zoo."""

from __future__ import annotations

from ..core import SequenceRegistry
from .binary_enots_wolley import BINARY_ENOTS_WOLLEY, BinaryEnotsWolleyGenerator
from .enots_wolley import ENOTS_WOLLEY, EnotsWolleyGenerator
from .forced_squarefree_enots_wolley import (
    FORCED_SQUAREFREE_ENOTS_WOLLEY,
    ForcedSquarefreeEnotsWolleyGenerator,
)
from .squarefree_semiprime_enots_wolley import (
    SQUAREFREE_SEMIPRIME_ENOTS_WOLLEY,
    SquarefreeSemiprimeEnotsWolleyGenerator,
)

__all__ = [
    "BINARY_ENOTS_WOLLEY",
    "BinaryEnotsWolleyGenerator",
    "ENOTS_WOLLEY",
    "EnotsWolleyGenerator",
    "FORCED_SQUAREFREE_ENOTS_WOLLEY",
    "ForcedSquarefreeEnotsWolleyGenerator",
    "SQUAREFREE_SEMIPRIME_ENOTS_WOLLEY",
    "SquarefreeSemiprimeEnotsWolleyGenerator",
    "register_builtins",
]


def register_builtins(registry: SequenceRegistry) -> None:
    registry.register(ENOTS_WOLLEY)
    registry.register(BINARY_ENOTS_WOLLEY)
    registry.register(FORCED_SQUAREFREE_ENOTS_WOLLEY)
    registry.register(SQUAREFREE_SEMIPRIME_ENOTS_WOLLEY)
