"""Compare retained-multiplier EW generation with the former generic filter path."""

from __future__ import annotations

import argparse
from time import perf_counter

from lex_earliest_seqs.zoo.every_kth_prime_only_enots_wolley import (
    EveryKthPrimeOnlyEnotsWolleyGenerator,
    EveryKthPrimeOnlyPolicy,
)
from lex_earliest_seqs.zoo.factor_restricted_enots_wolley import (
    FactorRestrictedEnotsWolleyGenerator,
)


def _elapsed(generator, count: int) -> float:
    start = perf_counter()
    generator.extend_to(count)
    return perf_counter() - start


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--count", type=int, default=10_000)
    parser.add_argument("--k", type=int, nargs="*", default=[2, 3, 4])
    args = parser.parse_args()

    if args.count < 2:
        parser.error("--count must be at least 2")

    print("k\tterms\tgeneric_s\tretained_s\tspeedup")
    for k in args.k:
        policy = EveryKthPrimeOnlyPolicy(k)

        generic = FactorRestrictedEnotsWolleyGenerator(policy=policy)
        generic_seconds = _elapsed(generic, args.count)

        retained = EveryKthPrimeOnlyEnotsWolleyGenerator(k=k)
        retained_seconds = _elapsed(retained, args.count)

        if generic.terms != retained.terms:
            raise RuntimeError(f"generator disagreement for k={k}")

        speedup = generic_seconds / retained_seconds
        print(
            f"{k}\t{args.count}\t{generic_seconds:.6f}\t"
            f"{retained_seconds:.6f}\t{speedup:.2f}x"
        )


if __name__ == "__main__":
    main()
