"""Built-in sequence zoo."""

from __future__ import annotations

from ..core import SequenceRegistry
from .binary_enots_wolley import BINARY_ENOTS_WOLLEY, BinaryEnotsWolleyGenerator
from .enots_wolley import ENOTS_WOLLEY, EnotsWolleyGenerator
from .every_kth_prime_enots_wolley import (
    EVERY_FOURTH_PRIME_ENOTS_WOLLEY,
    EVERY_KTH_PRIME_ENOTS_WOLLEY_DEFINITIONS,
    EVERY_SECOND_PRIME_ENOTS_WOLLEY,
    EVERY_THIRD_PRIME_ENOTS_WOLLEY,
    EveryKthPrimeEnotsWolleyGenerator,
    EveryKthPrimePolicy,
    ReferenceEveryKthPrimeEnotsWolleyGenerator,
    is_every_kth_prime,
    make_every_kth_prime_enots_wolley_definition,
    nth_prime,
    prime_index,
)
from .every_kth_prime_only_enots_wolley import (
    EVERY_FOURTH_PRIME_ONLY_ENOTS_WOLLEY,
    EVERY_KTH_PRIME_ONLY_ENOTS_WOLLEY_DEFINITIONS,
    EVERY_SECOND_PRIME_ONLY_ENOTS_WOLLEY,
    EVERY_THIRD_PRIME_ONLY_ENOTS_WOLLEY,
    EveryKthPrimeOnlyEnotsWolleyGenerator,
    EveryKthPrimeOnlyPolicy,
    ReferenceEveryKthPrimeOnlyEnotsWolleyGenerator,
    is_retained_prime,
    make_every_kth_prime_only_enots_wolley_definition,
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
from .sparse_prime_index_only_enots_wolley import (
    POWER_OF_TWO_INDEX_PRIME_ONLY_ENOTS_WOLLEY,
    SELF_POWER_INDEX_PRIME_ONLY_ENOTS_WOLLEY,
    SPARSE_PRIME_INDEX_ONLY_ENOTS_WOLLEY_DEFINITIONS,
    SQUARE_INDEX_PRIME_ONLY_ENOTS_WOLLEY,
    ReferenceSparsePrimeIndexOnlyEnotsWolleyGenerator,
    SparsePrimeIndexOnlyEnotsWolleyGenerator,
    SparsePrimeIndexOnlyPolicy,
    is_retained_prime_index,
    make_sparse_prime_index_only_enots_wolley_definition,
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
    "EVERY_FOURTH_PRIME_ENOTS_WOLLEY",
    "EVERY_FOURTH_PRIME_ONLY_ENOTS_WOLLEY",
    "EVERY_KTH_PRIME_ENOTS_WOLLEY_DEFINITIONS",
    "EVERY_KTH_PRIME_ONLY_ENOTS_WOLLEY_DEFINITIONS",
    "EVERY_SECOND_PRIME_ENOTS_WOLLEY",
    "EVERY_SECOND_PRIME_ONLY_ENOTS_WOLLEY",
    "EVERY_THIRD_PRIME_ENOTS_WOLLEY",
    "EVERY_THIRD_PRIME_ONLY_ENOTS_WOLLEY",
    "EWFactorPolicy",
    "EXACT_TRIPRIMARY_ENOTS_WOLLEY",
    "EnotsWolleyGenerator",
    "EveryKthPrimeEnotsWolleyGenerator",
    "EveryKthPrimeOnlyEnotsWolleyGenerator",
    "EveryKthPrimeOnlyPolicy",
    "EveryKthPrimePolicy",
    "FORCED_SQUAREFREE_ENOTS_WOLLEY",
    "FactorRestrictedEnotsWolleyDefinition",
    "FactorRestrictedEnotsWolleyGenerator",
    "ForcedSquarefreeEnotsWolleyGenerator",
    "POWER_OF_TWO_INDEX_PRIME_ONLY_ENOTS_WOLLEY",
    "PRIMARY_ENOTS_WOLLEY_DEFINITIONS",
    "ReferenceEveryKthPrimeEnotsWolleyGenerator",
    "ReferenceEveryKthPrimeOnlyEnotsWolleyGenerator",
    "ReferenceFactorRestrictedEnotsWolleyGenerator",
    "ReferenceSparsePrimeIndexOnlyEnotsWolleyGenerator",
    "SELF_POWER_INDEX_PRIME_ONLY_ENOTS_WOLLEY",
    "SPARSE_PRIME_INDEX_ONLY_ENOTS_WOLLEY_DEFINITIONS",
    "SQUAREFREE_EXACT_TRIPRIMARY_ENOTS_WOLLEY",
    "SQUAREFREE_SEMIPRIME_ENOTS_WOLLEY",
    "SQUAREFREE_TRIPRIMARY_ENOTS_WOLLEY",
    "SQUARE_INDEX_PRIME_ONLY_ENOTS_WOLLEY",
    "SparsePrimeIndexOnlyEnotsWolleyGenerator",
    "SparsePrimeIndexOnlyPolicy",
    "SquarefreeSemiprimeEnotsWolleyGenerator",
    "TRIPRIMARY_ENOTS_WOLLEY",
    "X000000_POLICY",
    "X000001_POLICY",
    "X000002_POLICY",
    "X000003_POLICY",
    "X000004_POLICY",
    "X000005_POLICY",
    "big_omega",
    "is_every_kth_prime",
    "is_retained_prime",
    "is_retained_prime_index",
    "make_every_kth_prime_enots_wolley_definition",
    "make_every_kth_prime_only_enots_wolley_definition",
    "make_factor_restricted_enots_wolley_definition",
    "make_sparse_prime_index_only_enots_wolley_definition",
    "nth_prime",
    "omega",
    "prime_index",
    "register_builtins",
]


def register_builtins(registry: SequenceRegistry) -> None:
    registry.register(ENOTS_WOLLEY)
    registry.register(BINARY_ENOTS_WOLLEY)
    registry.register(FORCED_SQUAREFREE_ENOTS_WOLLEY)
    for definition in PRIMARY_ENOTS_WOLLEY_DEFINITIONS:
        registry.register(definition)
    for definition in EVERY_KTH_PRIME_ENOTS_WOLLEY_DEFINITIONS:
        registry.register(definition)
    for definition in EVERY_KTH_PRIME_ONLY_ENOTS_WOLLEY_DEFINITIONS:
        registry.register(definition)
    for definition in SPARSE_PRIME_INDEX_ONLY_ENOTS_WOLLEY_DEFINITIONS:
        registry.register(definition)
