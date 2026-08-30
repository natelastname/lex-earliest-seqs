# Sequence identifiers

Sequences that have an OEIS entry use their six-digit `A`-number as the canonical ID, for example `A336957`.

Sequences that are not on the OEIS use a local six-digit `X`-number as the canonical ID, for example `X000000`. The `oeis` metadata field remains `None` until the sequence receives an actual OEIS entry.

## Local X-number assignments

| ID | Factor restriction |
| --- | --- |
| `X000000` | `omega(n) = 2`, squarefree |
| `X000001` | `omega(n) = 2`, multiplicity allowed |
| `X000002` | `omega(n) in {2, 3}`, multiplicity allowed |
| `X000003` | `omega(n) = 3`, multiplicity allowed |
| `X000004` | `omega(n) in {2, 3}`, squarefree |
| `X000005` | `omega(n) = 3`, squarefree |

Here `omega(n)` denotes the number of **distinct** prime factors. The independent `squarefree` flag controls whether repeated prime exponents are forbidden.
