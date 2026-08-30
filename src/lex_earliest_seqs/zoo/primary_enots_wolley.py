"""Finite-support-size Enots--Wolley family members with local X IDs."""

from __future__ import annotations

from .factor_restricted_enots_wolley import (
    EWFactorPolicy,
    make_factor_restricted_enots_wolley_definition,
)
from .squarefree_semiprime_enots_wolley import SQUAREFREE_SEMIPRIME_ENOTS_WOLLEY

# X000000 is the squarefree omega={2} member and has its own specialized
# closed-form generator. The remaining five requested members use the generic
# policy-driven generator for now.

X000001_POLICY = EWFactorPolicy(
    allowed_omega=frozenset({2}),
    squarefree=False,
)

X000002_POLICY = EWFactorPolicy(
    allowed_omega=frozenset({2, 3}),
    squarefree=False,
)

X000003_POLICY = EWFactorPolicy(
    allowed_omega=frozenset({2, 3}),
    squarefree=True,
)

X000004_POLICY = EWFactorPolicy(
    allowed_omega=frozenset({3}),
    squarefree=False,
)

X000005_POLICY = EWFactorPolicy(
    allowed_omega=frozenset({3}),
    squarefree=True,
)


BIPRIMARY_ENOTS_WOLLEY = make_factor_restricted_enots_wolley_definition(
    id="X000001",
    oeis=None,
    name="Biprimary Enots--Wolley",
    policy=X000001_POLICY,
    aliases=(
        "biprimary-ew",
        "biprimary-enots-wolley",
        "omega-2-ew",
    ),
    description=(
        "Lexicographically earliest sequence starting 1, 2 and obeying the "
        "Enots--Wolley rule while requiring every later term to have exactly "
        "two distinct prime factors; arbitrary positive prime multiplicities "
        "are allowed."
    ),
)

TRIPRIMARY_ENOTS_WOLLEY = make_factor_restricted_enots_wolley_definition(
    id="X000002",
    oeis=None,
    name="Triprimary Enots--Wolley",
    policy=X000002_POLICY,
    aliases=(
        "triprimary-ew",
        "triprimary-enots-wolley",
        "omega-2-3-ew",
    ),
    description=(
        "Lexicographically earliest sequence starting 1, 2 and obeying the "
        "Enots--Wolley rule while requiring every later term to have either "
        "two or three distinct prime factors; arbitrary positive prime "
        "multiplicities are allowed."
    ),
)

SQUAREFREE_TRIPRIMARY_ENOTS_WOLLEY = make_factor_restricted_enots_wolley_definition(
    id="X000003",
    oeis=None,
    name="Squarefree triprimary Enots--Wolley",
    policy=X000003_POLICY,
    aliases=(
        "squarefree-triprimary-ew",
        "squarefree-triprimary-enots-wolley",
        "squarefree-omega-2-3-ew",
    ),
    description=(
        "Lexicographically earliest sequence starting 1, 2 and obeying the "
        "Enots--Wolley rule while requiring every later term to be squarefree "
        "with either two or three distinct prime factors."
    ),
)

EXACT_TRIPRIMARY_ENOTS_WOLLEY = make_factor_restricted_enots_wolley_definition(
    id="X000004",
    oeis=None,
    name="Exact-triprimary Enots--Wolley",
    policy=X000004_POLICY,
    aliases=(
        "exact-triprimary-ew",
        "exact-triprimary-enots-wolley",
        "omega-3-ew",
    ),
    description=(
        "Lexicographically earliest sequence starting 1, 2 and obeying the "
        "Enots--Wolley rule while requiring every later term to have exactly "
        "three distinct prime factors; arbitrary positive prime multiplicities "
        "are allowed."
    ),
)

SQUAREFREE_EXACT_TRIPRIMARY_ENOTS_WOLLEY = (
    make_factor_restricted_enots_wolley_definition(
        id="X000005",
        oeis=None,
        name="Squarefree exact-triprimary Enots--Wolley",
        policy=X000005_POLICY,
        aliases=(
            "squarefree-exact-triprimary-ew",
            "squarefree-exact-triprimary-enots-wolley",
            "squarefree-omega-3-ew",
        ),
        description=(
            "Lexicographically earliest sequence starting 1, 2 and obeying the "
            "Enots--Wolley rule while requiring every later term to be "
            "squarefree with exactly three distinct prime factors."
        ),
    )
)

PRIMARY_ENOTS_WOLLEY_DEFINITIONS = (
    SQUAREFREE_SEMIPRIME_ENOTS_WOLLEY,
    BIPRIMARY_ENOTS_WOLLEY,
    TRIPRIMARY_ENOTS_WOLLEY,
    SQUAREFREE_TRIPRIMARY_ENOTS_WOLLEY,
    EXACT_TRIPRIMARY_ENOTS_WOLLEY,
    SQUAREFREE_EXACT_TRIPRIMARY_ENOTS_WOLLEY,
)
