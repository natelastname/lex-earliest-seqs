# Sequence identifiers

Sequences that have an OEIS entry use their six-digit `A`-number as the canonical ID, for example `A336957`.

Sequences that are not on the OEIS use a local six-digit `X`-number as the canonical ID, for example `X000000`. The `oeis` metadata field remains `None` until the sequence receives an actual OEIS entry.

## Local X-number assignments

| ID | Factor restriction |
| --- | --- |
| `X000000` | `omega(n) = 2`, squarefree |
| `X000001` | `omega(n) = 2`, multiplicity allowed |
| `X000002` | `omega(n) in {2, 3}`, multiplicity allowed |
| `X000003` | `omega(n) in {2, 3}`, squarefree |
| `X000004` | `omega(n) = 3`, multiplicity allowed |
| `X000005` | `omega(n) = 3`, squarefree |
| `X000006` | divisible by at least one prime `p_j` with `2 | j` |
| `X000007` | divisible by at least one prime `p_j` with `3 | j` |
| `X000008` | divisible by at least one prime `p_j` with `4 | j` |

Here `omega(n)` denotes the number of **distinct** prime factors. The independent `squarefree` flag controls whether repeated prime exponents are forbidden.

The `X000006`--`X000008` family uses the general every-`k`-th-prime rule: for a parameter `k`, the distinguished prime coordinates are `p_k, p_{2k}, p_{3k}, ...`, and a term is eligible exactly when at least one of those primes divides it. All other primes remain legal as cofactors; they simply do not make a term eligible by themselves. Only `k = 2, 3, 4` are currently registered, although the generator accepts arbitrary positive `k`.
