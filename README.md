# lex-earliest-seqs

A small research library for defining, computing, caching, and inspecting
lexicographically earliest sequences.

The package deliberately distinguishes a **sequence position** from the
**ambient rank of an object**. A sequence chooses objects from an ordered
object space; the 17th sequence term need not be the 17th admissible object.

## Design

A sequence definition supplies:

- a stable ID and name (optionally an OEIS number),
- an ordered ambient object space,
- a factory for a mutable stateful generator,
- a generator/definition version for cache compatibility, and
- zero or more incidence projections used only for chronology tables.

A generator has only one required contract:

```python
class MyGenerator:
    terms: list[MyObject]

    def extend_to(self, count: int) -> None:
        ...
```

The **entire generator object is pickled**. Sequence-specific generators are
therefore free to retain heaps, bitsets, factor tables, occurrence indexes, or
other expensive computation state without teaching the generic package how to
serialize them. Python arbitrary-precision integers are preserved directly.
Caches are trusted local research artifacts; pickle security is intentionally
not treated as a package concern.

Incidence projections are separate from generation. Changing a support
function or chronology-table renderer cannot change the sequence.

## Example

```python
from lex_earliest_seqs import open_run, registry
from lex_earliest_seqs.incidence import build_incidence_table, render_text

sequence = registry.resolve("ew")
run = open_run(sequence)
run.ensure(1000)

projection = sequence.projections["prime-exponents"]
table = build_incidence_table(
    run.records(0, 30),
    projection=projection,
)
print(render_text(table))
```

The initial zoo contains Enots--Wolley (OEIS A336957) and its binary-support
analogue (OEIS A338833). Both generators are intentionally transparent rather
than optimized; they serve as reference implementations and demonstrations of
pickleable stateful generators with different incidence projections.

## Defining a sequence

```python
from dataclasses import dataclass, field

from lex_earliest_seqs import PositiveIntegers, SequenceDefinition
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
    generator_factory=Generator,
    generator_version=1,
    object_space=PositiveIntegers(),
    projections={"bits": binary_digit_projection()},
)
```

A boolean support function can be adapted directly with
`boolean_support_projection`; weighted supports can return sparse coordinate
mappings through `IncidenceProjection`.

## CLI

```console
lex-earliest-seqs list
lex-earliest-seqs info ew
lex-earliest-seqs compute ew 10000
lex-earliest-seqs terms ew 20
lex-earliest-seqs table ew 30 --projection prime-exponents
lex-earliest-seqs table A338833 30 --projection binary-digits
lex-earliest-seqs table ew 30 --format markdown
```

By default pickles are stored under
`$XDG_CACHE_HOME/lex-earliest-seqs` or `~/.cache/lex-earliest-seqs`.
Writes use a temporary file followed by an atomic replace.

## Development

```console
uv sync
uv run pytest
```

## License

MIT / Expat
