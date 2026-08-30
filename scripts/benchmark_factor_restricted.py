"""Benchmark finite-omega EW generators against correctness/performance baselines.

Examples:

    uv run python scripts/benchmark_factor_restricted.py --terms 100 500 1000

    uv run python scripts/benchmark_factor_restricted.py \
      --ids X000001 X000002 X000004 \
      --terms 1000 5000 \
      --compare-retirement

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
    FactorRestrictedEnotsWolleyGenerator,
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


class LazyRetirementFactorRestrictedEnotsWolleyGenerator(
    FactorRestrictedEnotsWolleyGenerator
):
    """A/B benchmark baseline reproducing pre-eager retirement behavior."""

    def _retire_used_value(self, value: int) -> int:
        return 0


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


def _assert_equal(left, right, *, label: str) -> None:
    if left.terms == right.terms:
        return
    for index, (left_term, right_term) in enumerate(
        zip(left.terms, right.terms, strict=True),
        start=1,
    ):
        if left_term != right_term:
            raise AssertionError(
                f"{label} diverges at term {index}: "
                f"left={left_term}, right={right_term}"
            )
    raise AssertionError(f"{label} length mismatch")


def _verified_factories(
    definition: FactorRestrictedEnotsWolleyDefinition,
    count: int,
):
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
    _assert_equal(optimized, reference, label=f"{definition.id} optimized/reference")

    return optimized_factory, reference_factory


def _verified_retirement_factories(
    definition: FactorRestrictedEnotsWolleyDefinition,
    count: int,
):
    policy = definition.factor_policy
    assert policy is not None

    eager_factory = lambda: FactorRestrictedEnotsWolleyGenerator(policy=policy)
    lazy_factory = lambda: LazyRetirementFactorRestrictedEnotsWolleyGenerator(
        policy=policy
    )

    _clear_factor_caches()
    eager = eager_factory()
    eager.extend_to(count)
    _clear_factor_caches()
    lazy = lazy_factory()
    lazy.extend_to(count)
    _assert_equal(eager, lazy, label=f"{definition.id} eager/lazy retirement")

    return eager_factory, lazy_factory


def _print_reference_benchmark(args) -> None:
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


def _print_retirement_benchmark(args) -> None:
    print()
    print("### Eager retirement A/B")
    print()
    print("| sequence | terms | eager (s) | lazy (s) | lazy/eager | verified |")
    print("| --- | ---: | ---: | ---: | ---: | :---: |")

    for sequence_id in args.ids:
        if sequence_id == "X000000":
            # X000000 uses its theorem-based closed-form generator in production;
            # this A/B comparison is specifically about the generic stream engine.
            continue

        definition = registry.resolve(sequence_id)
        if not isinstance(definition, FactorRestrictedEnotsWolleyDefinition):
            raise TypeError(f"{sequence_id} is not a factor-restricted EW definition")

        for count in args.terms:
            eager_factory, lazy_factory = _verified_retirement_factories(
                definition,
                count,
            )
            eager_seconds = _time_fresh(eager_factory, count, args.repeat)
            lazy_seconds = _time_fresh(lazy_factory, count, args.repeat)
            ratio = lazy_seconds / eager_seconds
            print(
                f"| {sequence_id} | {count:,} | {eager_seconds:.6f} | "
                f"{lazy_seconds:.6f} | {ratio:.2f}x | yes |"
            )


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
    parser.add_argument(
        "--compare-retirement",
        action="store_true",
        help="also compare eager cross-stream retirement with the prior lazy strategy",
    )
    args = parser.parse_args()

    if args.repeat < 1:
        raise ValueError("--repeat must be positive")
    if any(count < 2 for count in args.terms):
        raise ValueError("all --terms values must be at least 2")

    _print_reference_benchmark(args)
    if args.compare_retirement:
        _print_retirement_benchmark(args)


if __name__ == "__main__":
    main()
