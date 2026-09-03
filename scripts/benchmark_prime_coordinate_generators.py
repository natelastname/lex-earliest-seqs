"""Benchmark the prime-coordinate EW generator optimization layers."""

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


class SingleTableBaselineGenerator(EveryKthPrimeOnlyEnotsWolleyGenerator):
    """Retained-monoid generator without the odd-table fast path."""

    def _multiplier_table_for_forbidden(self, forbidden_radical):
        return self.multiplier_values, self.multiplier_successors


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

    print(
        "k\tterms\tgeneric_s\tsingle_table_s\todd_fast_s\t"
        "vs_generic\tvs_single_table"
    )
    for k in args.k:
        policy = EveryKthPrimeOnlyPolicy(k)

        generic = FactorRestrictedEnotsWolleyGenerator(policy=policy)
        generic_seconds = _elapsed(generic, args.count)

        single_table = SingleTableBaselineGenerator(k=k)
        single_table_seconds = _elapsed(single_table, args.count)

        optimized = EveryKthPrimeOnlyEnotsWolleyGenerator(k=k)
        optimized_seconds = _elapsed(optimized, args.count)

        if generic.terms != single_table.terms or generic.terms != optimized.terms:
            raise RuntimeError(f"generator disagreement for k={k}")

        print(
            f"{k}\t{args.count}\t{generic_seconds:.6f}\t"
            f"{single_table_seconds:.6f}\t{optimized_seconds:.6f}\t"
            f"{generic_seconds / optimized_seconds:.2f}x\t"
            f"{single_table_seconds / optimized_seconds:.2f}x"
        )


if __name__ == "__main__":
    main()
