from lex_earliest_seqs import registry


def test_generic_factor_family_cache_versions_include_eager_retirement():
    # X000000 has its own theorem-based generator and is unaffected by the
    # persistent-stream retirement change. The five generic stream generators
    # moved from version 2 to 3 so old caches cannot silently retain incomplete
    # cross-stream retirement history.
    assert registry.resolve("X000000").generator_version == 2
    for sequence_id in (
        "X000001",
        "X000002",
        "X000003",
        "X000004",
        "X000005",
    ):
        assert registry.resolve(sequence_id).generator_version == 3
