# Benchmarks

## Factor-restricted Enots--Wolley family

The production generator for the finite-omega family uses persistent, path-compressed candidate streams. The old direct positive-integer scanner is retained as `ReferenceFactorRestrictedEnotsWolleyGenerator` and serves as both a correctness oracle and a benchmark baseline.

Run the benchmark with:

```console
uv run python scripts/benchmark_factor_restricted.py --terms 100 500 1000
```

For a longer biprimary-only run:

```console
uv run python scripts/benchmark_factor_restricted.py \
  --ids X000001 \
  --terms 1000 5000 10000 \
  --repeat 3
```

The benchmark:

1. constructs fresh generators rather than loading sequence caches;
2. clears the shared factorization/support caches between timed runs;
3. verifies exact optimized/reference prefix equality before timing each cell;
4. reports the median of the requested number of repetitions;
5. reports `reference / optimized` as the speedup factor.

`X000000` is included by default even though its production generator is the stronger theorem-based closed-form implementation. `X000001` through `X000005` exercise the generic persistent-stream engine.

Benchmark results are intentionally not used as correctness tests or hard performance gates: absolute timings are machine-dependent, while the optimized/reference equality checks belong in the pytest suite.
