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
from .full_return_enots_wolley import (
    FULL_RETURN_ENOTS_WOLLEY_DEFINITIONS,
    FULL_RETURN_EW_2_3,
    FULL_RETURN_EW_2_5,
    FULL_RETURN_EW_3_5,
    FullReturnEnotsWolleyGenerator,
    ReferenceFullReturnEnotsWolleyGenerator,
    full_return_candidate_allowed,
    full_return_restriction_active,
    make_full_return_enots_wolley_definition,
    target_free,
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
    "FULL_RETURN_ENOTS_WOLLEY_DEFINITIONS",
    "FULL_RETURN_EW_2_3",
    "FULL_RETURN_EW_2_5",
    "FULL_RETURN_EW_3_5",
    "FactorRestrictedEnotsWolleyDefinition",
    "FactorRestrictedEnotsWolleyGenerator",
    "ForcedSquarefreeEnotsWolleyGenerator",
    "FullReturnEnotsWolleyGenerator",
    "PRIMARY_ENOTS_WOLLEY_DEFINITIONS",
    "ReferenceEveryKthPrimeEnotsWolleyGenerator",
    "ReferenceEveryKthPrimeOnlyEnotsWolleyGenerator",
    "ReferenceFactorRestrictedEnotsWolleyGenerator",
    "ReferenceFullReturnEnotsWolleyGenerator",
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
    "full_return_candidate_allowed",
    "full_return_restriction_active",
    "is_every_kth_prime",
    "is_retained_prime",
    "make_every_kth_prime_enots_wolley_definition",
    "make_every_kth_prime_only_enots_wolley_definition",
    "make_factor_restricted_enots_wolley_definition",
    "make_full_return_enots_wolley_definition",
    "nth_prime",
    "omega",
    "prime_index",
    "register_builtins",
    "target_free",
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
    for definition in FULL_RETURN_ENOTS_WOLLEY_DEFINITIONS:
        registry.register(definition)
