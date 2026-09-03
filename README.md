# lex-earliest-seqs

Research tools for defining, computing, caching, and inspecting
**lexicographically earliest sequences**.

The package is designed for sequences whose next term is chosen greedily from
an ordered ambient object space. It keeps sequence positions separate from
ambient-object ranks, supports stateful high-performance generators, persists
entire generator states with pickle, and provides generic incidence chronology
tables for supports such as prime factors or binary digits.

Python 3.12+ is required.
