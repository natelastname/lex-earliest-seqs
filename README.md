# lex-earliest-seqs

Research tools for defining, computing, caching, and inspecting
**lexicographically earliest sequences**.

The package is designed for sequences whose next term is chosen greedily from
an ordered ambient object space. It keeps sequence positions separate from
ambient-object ranks, supports stateful high-performance generators, persists
entire generator states with pickle, and provides generic incidence chronology
tables for supports such as prime factors or binary digits.

Python 3.12+ is required.

## Major features

### Sequence position vs. ambient object

The package deliberately distinguishes:

- **position**: zero-based position in the generated sequence;
- **subscript**: mathematical sequence index (`position + offset`);
- **object rank**: zero-based rank in the ambient admissible object space;
- **value**: the actual object selected by the sequence.

For example, the 17th sequence term need not be the 17th admissible integer.
This distinction lets the same framework work with positive integers,
nonnegative integers, or other ordered countable object spaces.

### Stateful, pickleable generators

A generator has a deliberately small contract:

```python
class MyGenerator:
    terms: list[MyObject]

    def extend_to(self, count: int) -> None:
        ...
```

The **entire generator object is persisted**. Sequence-specific implementations
can therefore retain whatever continuation state makes computation fast:
heaps, occurrence indexes, candidate frontiers, factor tables, or other data
structures. Derived structures can be omitted from the pickle and rebuilt
lazily after loading.

Caches are trusted local research artifacts. Pickle security against untrusted
files is intentionally not a package concern.

### Persistent computation cache

`open_run()` automatically loads and saves generator pickles. By default they
live under:

```text
$XDG_CACHE_HOME/lex-earliest-seqs/
```

or, when `XDG_CACHE_HOME` is unset:

```text
~/.cache/lex-earliest-seqs/
```

Writes use a temporary file followed by an atomic replace. Sequence definitions
carry generator and definition versions so incompatible caches are rejected
rather than silently reused.

### Incidence chronology tables

Generation and incidence/display logic are separate. A sequence can expose one
or more named **incidence projections**, such as:

- prime exponents / prime support;
- binary digit support;
- a custom Boolean support function;
- a custom sparse weighted support.

Chronology tables can be rendered as text, Markdown, JSON, CSV, or TSV. Changing
a projection or renderer cannot change sequence generation.

### Built-in sequence zoo

| ID | Name | Useful aliases | Projection |
| --- | --- | --- | --- |
| `A336957` | Enots--Wolley | `ew`, `enots-wolley` | `prime-exponents` |
| `A338833` | Binary Enots--Wolley | `bew`, `binary-ew` | `binary-digits` |
| `A399457` | Forced-squarefree Enots--Wolley | `squarefree-ew`, `forced-squarefree-ew` | `prime-exponents` |

The built-ins use sequence-specific candidate enumeration rather than relying
on a generic brute-force scanner:

- **A336957** partitions locally eligible integers into disjoint prime-indexed
  streams, merges them with a heap, and persists lazy per-prime multiplier
  successor maps so recurring streams jump over products already used earlier;
- **A338833** uses an exact bit successor to jump directly between locally
  admissible binary candidates;
- **A399457** merges squarefree candidate streams and persists per-stream
  cofactor frontiers so late computation skips historically exhausted prefixes.

## Installation / development

```console
uv sync
```

Run the CLI through the project environment with:

```console
uv run lex-earliest-seqs list
```

The installed console command is `lex-earliest-seqs`.

## Command-line interface

Sequence arguments accept the canonical ID, OEIS number, registered name, or
alias.

### List the built-in sequences

```console
lex-earliest-seqs list
```

### Inspect sequence metadata

```console
lex-earliest-seqs info ew
lex-earliest-seqs info A399457
```

This prints the ID, name, OEIS number, offset, object space, cache versions,
registered projections, cached term count, and description. If the sequence has
no cache, `cached terms` is reported as `0`. Use `--cache-dir` to inspect a
non-default cache directory:

```console
lex-earliest-seqs info ew --cache-dir ./sequence-cache
```

### Compute/cache a prefix

```console
lex-earliest-seqs compute ew 100000
lex-earliest-seqs compute A399457 1000000
```

`compute` ensures that at least the requested number of terms are present in
the generator cache.

Useful cache controls are:

```console
# Ignore any existing cache and start from a fresh generator.
lex-earliest-seqs compute ew 10000 --refresh

# Do not load or save a cache.
lex-earliest-seqs compute ew 10000 --no-cache

# Store caches in another directory.
lex-earliest-seqs compute ew 10000 --cache-dir ./sequence-cache

# Suppress progress messages.
lex-earliest-seqs compute ew 10000 --no-progress
```

Cache-loading and computation progress are written to **stderr**, so sequence
or table data on stdout remains clean for piping.

### Print or export terms

```console
lex-earliest-seqs terms ew 20
lex-earliest-seqs terms A399457 20
```

Without an output file, the output contains mathematical subscript and value,
separated by a tab. Slices use a zero-based sequence position:

```console
lex-earliest-seqs terms ew 25 --start-position 1000
```

Terms can also be exported as CSV or Parquet:

```console
lex-earliest-seqs terms ew 1000000 -o ew.parquet --format parquet
lex-earliest-seqs terms ew 1000000 -o ew.csv --format csv
```

`--output` is equivalent to `-o`. When `--format` is omitted, `.csv`,
`.parquet`, and `.pq` filename suffixes select the corresponding format. Both
file formats use stable `subscript` and `value` columns. Parquet output is
written in bounded batches with Zstd compression, so the export step does not
construct a second full sequence-sized record list in memory.

For example, an integer-valued sequence export can be passed directly to
`massive-scatter`:

```console
massive-scatter build ew.parquet ew.msplot --x subscript --y value
```

The same `--cache-dir`, `--refresh`, `--no-cache`, and `--no-progress` controls
available to `compute` also work for `terms`.

### Print an incidence chronology table

```console
lex-earliest-seqs table ew 30
lex-earliest-seqs table A338833 30 --projection binary-digits
lex-earliest-seqs table A399457 50 --projection prime-exponents
```

When a sequence has exactly one projection, `--projection` may be omitted.

Available output formats are:

```console
lex-earliest-seqs table ew 30 --format text
lex-earliest-seqs table ew 30 --format markdown
lex-earliest-seqs table ew 30 --format json
lex-earliest-seqs table ew 30 --format csv
lex-earliest-seqs table ew 30 --format tsv
```

Other useful table options:

```console
# Start at zero-based sequence position 1000.
lex-earliest-seqs table ew 30 --start-position 1000

# Include every feature column through the largest encountered feature,
# rather than only columns actually used in the selected rows.
lex-earliest-seqs table ew 30 --columns through-largest

# Control panel splitting for text output.
lex-earliest-seqs table ew 50 --width 160
```

`--columns` accepts `used` (the default) or `through-largest`.

## Python API

The main public API is exported directly from `lex_earliest_seqs`.

### Resolve and compute a sequence

```python
from lex_earliest_seqs import open_run, registry

definition = registry.resolve("A399457")
run = open_run(definition)

run.ensure(100_000)

print(run.terms[:20])
print(run.terms[99_999])
```
