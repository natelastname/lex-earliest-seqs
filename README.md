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

Aliases and OEIS IDs resolve to the same immutable `SequenceDefinition`:

```python
ew = registry.resolve("ew")
binary = registry.resolve("bew")
squarefree = registry.resolve("squarefree-ew")
```

Iterate over the registry to discover definitions programmatically:

```python
from lex_earliest_seqs import registry

for definition in registry:
    print(definition.id, definition.name, tuple(definition.projections))
```

### Positions, subscripts, and records

`SequenceRun` provides explicit accessors for the different coordinate systems:

```python
from lex_earliest_seqs import open_run, registry

run = open_run(registry.resolve("ew"))

# Zero-based sequence position.
value = run.at_position(99)

# Mathematical sequence subscript, respecting the definition's offset.
value = run.at_subscript(100)

record = run.record_at_position(99)
print(record.position)
print(record.subscript)
print(record.object_rank)
print(record.value)
```

Get a range of full `TermRecord` objects with:

```python
records = run.records(100, 130)
```

### Cache controls from Python

Default cached run:

```python
run = open_run(registry.resolve("ew"))
```

Disable caching entirely:

```python
run = open_run(registry.resolve("ew"), use_cache=False)
```

Use another cache directory:

```python
from pathlib import Path

run = open_run(
    registry.resolve("ew"),
    cache_dir=Path("./sequence-cache"),
)
```

Use one exact cache file:

```python
run = open_run(
    registry.resolve("ew"),
    cache_path="./ew.pkl",
)
```

Ignore an existing cache and construct a fresh generator:

```python
run = open_run(registry.resolve("ew"), refresh=True)
```

`cache_dir` and `cache_path` are mutually exclusive.

### Computation progress from Python

`SequenceRun.ensure()` accepts a callback receiving `(current, target)`:

```python
from lex_earliest_seqs import open_run, registry

run = open_run(registry.resolve("A399457"))


def progress(current: int, target: int) -> None:
    print(f"{current:,}/{target:,}")


run.ensure(1_000_000, progress=progress)
```

Progress-aware computation extends opaque generators in batches. You can
control the batch size explicitly with `progress_chunk_size`.

Cache loading has a separate byte-progress callback:

```python
run = open_run(
    registry.resolve("ew"),
    load_progress=lambda current, total: print(current, total),
)
```

### Build chronology tables from Python

```python
from lex_earliest_seqs import (
    build_incidence_table,
    open_run,
    registry,
    render_text,
)

sequence = registry.resolve("ew")
run = open_run(sequence)
run.ensure(30)

projection = sequence.projections["prime-exponents"]
table = build_incidence_table(
    run.records(0, 30),
    projection=projection,
)

print(render_text(table))
```

Other public renderers include:

```python
from lex_earliest_seqs import (
    render_delimited,
    render_json,
    render_markdown,
)

markdown = render_markdown(table)
json_text = render_json(table)
csv_text = render_delimited(table, delimiter=",")
tsv_text = render_delimited(table, delimiter="\t")
```

### Define and register a new sequence

A new sequence supplies an ordered object space, a stateful generator factory,
metadata, and optionally incidence projections.

```python
from dataclasses import dataclass, field

from lex_earliest_seqs import (
    PositiveIntegers,
    SequenceDefinition,
    registry,
)
from lex_earliest_seqs.projections import binary_digit_projection


@dataclass
class Generator:
    terms: list[int] = field(default_factory=lambda: [1])

    def extend_to(self, count: int) -> None:
        while len(self.terms) < count:
            self.terms.append(self.terms[-1] + 1)


MY_SEQUENCE = SequenceDefinition(
    id="my-sequence",
    name="My sequence",
    aliases=("mine",),
    generator_factory=Generator,
    generator_version=1,
    definition_version=1,
    offset=1,
    object_space=PositiveIntegers(),
    projections={"bits": binary_digit_projection()},
)

registry.register(MY_SEQUENCE)
```

The generator may contain arbitrary pickleable continuation state. If a code
change makes previously persisted generator state incompatible, increment
`generator_version`. If the mathematical/metadata definition changes in a way
that invalidates caches, increment `definition_version`.

Boolean support functions can be adapted with
`boolean_support_projection`; weighted supports can be represented as sparse
coordinate mappings through `IncidenceProjection`.

## One-off migration of the old EW cache

The temporary migration helper converts the historical `enots-wolley-2`
A336957 term-list pickle into the native stateful-generator cache **without
recomputing or replaying the sequence prefix**:

```console
uv run python scripts/migrate_legacy_ew_cache.py
```

By default it reads:

```text
~/.cache/enots-wolley-2/terms-v1.pkl
```

and writes:

```text
~/.cache/lex-earliest-seqs/A336957.pkl
```

If the native target already exists, pass `--force`. This migration code is
intentionally temporary and can be removed after the research cache has been
converted and verified.

## Development

```console
uv sync
uv run pytest
```

## License

MIT / Expat
