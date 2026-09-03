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
| `X000006` | contains at least one prime `p_j` with `2 | j`; other prime factors allowed |
| `X000007` | contains at least one prime `p_j` with `3 | j`; other prime factors allowed |
| `X000008` | contains at least one prime `p_j` with `4 | j`; other prime factors allowed |
| `X000009` | every prime factor lies in `p_1, p_3, p_5, ...` (`k = 2`) |
| `X000010` | every prime factor lies in `p_1, p_4, p_7, ...` (`k = 3`) |
| `X000011` | every prime factor lies in `p_1, p_5, p_9, ...` (`k = 4`) |

Here `omega(n)` denotes the number of **distinct** prime factors. The independent `squarefree` flag controls whether repeated prime exponents are forbidden.

`X000006`--`X000008` are the earlier **contains-a-distinguished-prime** family. For parameter `k`, the distinguished prime coordinates are `p_k, p_{2k}, p_{3k}, ...`; a term is eligible when at least one such prime divides it, while all other primes remain legal as cofactors.

`X000009`--`X000011` are the **prime-coordinate-restricted** family. For parameter `k`, only coordinates `p_1, p_{1+k}, p_{1+2k}, ...` exist. Every prime factor of every generated noninitial term must come from this retained subsequence. Thus for `k = 3`, the allowed primes begin `2, 7, 17, 29, ...`, so no term is divisible by `3` or `5`. Both family generators accept arbitrary positive `k`; only `k = 2, 3, 4` are registered.
