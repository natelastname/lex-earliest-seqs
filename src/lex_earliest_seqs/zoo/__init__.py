"""Built-in sequence zoo."""

from __future__ import annotations

from ..core import SequenceRegistry
from .binary_enots_wolley import BINARY_ENOTS_WOLLEY, BinaryEnotsWolleyGenerator
from .enots_wolley import ENOTS_WOLLEY, EnotsWolleyGenerator
from .even_index_prime_enots_wolley import (
    EVEN_INDEX_PRIME_ENOTS_WOLLEY,
    EvenIndexPrimeEnotsWolleyGenerator,
    EvenIndexPrimePolicy,
    ReferenceEvenIndexPrimeEnotsWolleyGenerator,
    is_even_index_prime,
)
from .factor_restricted_enots_wolley import (
    EWFactorPolicy,
    FactorRestrictedEnotsWolleyDefinition,
    FactorRestrictedEnotsWolleyGenerator,
    ReferenceFactorRestrictedEnotsWolleyGenerator,
    big_omega,
    make_factor_restricted_enots_wolley_definition,
    omega,
)
from .forced_squarefree_enots_wolley import (
    FORCED_SQUAREFREE_ENOTS_WOLLEY,
    ForcedSquarefreeEnotsWolleyGenerator,
)
from .primary_enots_wolley import (
    BIPRIMARY_ENOTS_WOLLEY,
    EXACT_TRIPRIMARY_ENOTS_WOLLEY,
    PRIMARY_ENOTS_WOLLEY_DEFINITIONS,
    SQUAREFREE_EXACT_TRIPRIMARY_ENOTS_WOLLEY,
    SQUAREFREE_TRIPRIMARY_ENOTS_WOLLEY,
    TRIPRIMARY_ENOTS_WOLLEY,
    X000001_POLICY,
    X000002_POLICY,
    X000003_POLICY,
    X000004_POLICY,
    X000005_POLICY,
)
from .squarefree_semiprime_enots_wolley import (
    SQUAREFREE_SEMIPRIME_ENOTS_WOLLEY,
    X000000_POLICY,
    SquarefreeSemiprimeEnotsWolleyGenerator,
)

__all__ = [
    "BINARY_ENOTS_WOLLEY",
    "BIPRIMARY_ENOTS_WOLLEY",
    "BinaryEnotsWolleyGenerator",
    "ENOTS_WOLLEY",
    "EVEN_INDEX_PRIME_ENOTS_WOLLEY",
    "EWFactorPolicy",
    "EXACT_TRIPRIMARY_ENOTS_WOLLEY",
    "EnotsWolleyGenerator",
    "EvenIndexPrimeEnotsWolleyGenerator",
    "EvenIndexPrimePolicy",
    "FORCED_SQUAREFREE_ENOTS_WOLLEY",
    "FactorRestrictedEnotsWolleyDefinition",
    "FactorRestrictedEnotsWolleyGenerator",
    "ForcedSquarefreeEnotsWolleyGenerator",
    "PRIMARY_ENOTS_WOLLEY_DEFINITIONS",
    "ReferenceEvenIndexPrimeEnotsWolleyGenerator",
    "ReferenceFactorRestrictedEnotsWolleyGenerator",
    "SQUAREFREE_EXACT_TRIPRIMARY_ENOTS_WOLLEY",
    "SQUAREFREE_SEMIPRIME_ENOTS_WOLLEY",
    "SQUAREFREE_TRIPRIMARY_ENOTS_WOLLEY",
    "SquarefreeSemiprimeEnotsWolleyGenerator",
    "TRIPRIMARY_ENOTS_WOLLEY",
    "X000000_POLICY",
    "X000001_POLICY",
    "X000002_POLICY",
    "X000003_POLICY",
    "X000004_POLICY",
    "X000005_POLICY",
    "big_omega",
    "is_even_index_prime",
    "make_factor_restricted_enots_wolley_definition",
    "omega",
    "register_builtins",
]


def register_builtins(registry: SequenceRegistry) -> None:
    registry.register(ENOTS_WOLLEY)
    registry.register(BINARY_ENOTS_WOLLEY)
    registry.register(FORCED_SQUAREFREE_ENOTS_WOLLEY)
    for definition in PRIMARY_ENOTS_WOLLEY_DEFINITIONS:
        registry.register(definition)
    registry.register(EVEN_INDEX_PRIME_ENOTS_WOLLEY)
