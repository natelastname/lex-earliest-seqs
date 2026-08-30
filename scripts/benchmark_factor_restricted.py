"""Benchmark optimized finite-omega EW generators against the direct oracle.

Example:

    uv run python scripts/benchmark_factor_restricted.py --terms 100 500 1000

The script always verifies exact sequence equality before reporting timing data.
"""

from __future__ import annotations

import argparse
import gc
import statistics
import time
from collections.abc import Callable

from lex_earliest_seqs import registry
from lex_earliest_seqs.projections import prime_factorization
from lex_earliest_seqs.zoo.enots_wolley import prime_support
from lex_earliest_seqs.zoo.factor_restricted_enots_wolley import (
    FactorRestrictedEnotsWolleyDefinition,
    ReferenceFactorRestrictedEnotsWolleyGenerator,
)

DEFAULT_IDS = (
    "X000000",
    "X000001",
    "X000002",
    "X000003",
    "X000004",
    "X000005",
)


def _clear_factor_caches() -> None:
    prime_factorization.cache_clear()
    prime_support.cache_clear()


def _time_fresh(factory: Callable[[], object], count: int, repeat: int) -> float:
    samples: list[float] = []
    for _ in range(repeat):
        _clear_factor_caches()
        gc.collect()
        generator = factory()
        start = time.perf_counter()
        generator.extend_to(count)
        samples.append(time.perf_counter() - start)
    return statistics.median(samples)


def _verified_factories(definition: FactorRestrictedEnotsWolleyDefinition, count: int):
    policy = definition.factor_policy
    assert policy is not None

    optimized_factory = definition.generator_factory
    reference_factory = lambda: ReferenceFactorRestrictedEnotsWolleyGenerator(
        policy=policy
    )

    _clear_factor_caches()
    optimized = optimized_factory()
    optimized.extend_to(count)
    _clear_factor_caches()
    reference = reference_factory()
    reference.extend_to(count)
    if optimized.terms != reference.terms:
        for index, (left, right) in enumerate(
            zip(optimized.terms, reference.terms, strict=True),
            start=1,
        ):
            if left != right:
                raise AssertionError(
                    f"{definition.id} diverges at term {index}: "
                    f"optimized={left}, reference={right}"
                )
        raise AssertionError(f"{definition.id} optimized/reference length mismatch")

    return optimized_factory, reference_factory


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--terms",
        nargs="+",
        type=int,
        default=[100, 500, 1_000],
        help="prefix sizes to benchmark (default: 100 500 1000)",
    )
    parser.add_argument(
        "--ids",
        nargs="+",
        default=list(DEFAULT_IDS),
        help="sequence IDs to benchmark",
    )
    parser.add_argument(
        "--repeat",
        type=int,
        default=3,
        help="timed fresh runs per cell; median is reported (default: 3)",
    )
    args = parser.parse_args()

    if args.repeat < 1:
        raise ValueError("--repeat must be positive")
    if any(count < 2 for count in args.terms):
        raise ValueError("all --terms values must be at least 2")

    print(
        "| sequence | terms | optimized (s) | reference (s) | speedup | verified |"
    )
    print("| --- | ---: | ---: | ---: | ---: | :---: |")

    for sequence_id in args.ids:
        definition = registry.resolve(sequence_id)
        if not isinstance(definition, FactorRestrictedEnotsWolleyDefinition):
            raise TypeError(f"{sequence_id} is not a factor-restricted EW definition")

        for count in args.terms:
            optimized_factory, reference_factory = _verified_factories(definition, count)
            optimized_seconds = _time_fresh(optimized_factory, count, args.repeat)
            reference_seconds = _time_fresh(reference_factory, count, args.repeat)
            speedup = reference_seconds / optimized_seconds
            print(
                f"| {sequence_id} | {count:,} | {optimized_seconds:.6f} | "
                f"{reference_seconds:.6f} | {speedup:.2f}x | yes |"
            )


if __name__ == "__main__":
    main()
