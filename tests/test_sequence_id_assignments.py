from lex_earliest_seqs import registry


def test_local_sequence_ids_are_contiguous_through_x000014():
    local_ids = [definition.id for definition in registry if definition.id.startswith("X")]
    assert local_ids == [f"X{index:06d}" for index in range(15)]


def test_repaired_sequence_families_have_canonical_ids_and_names():
    expected = {
        "X000006": "Every-second-prime Enots--Wolley",
        "X000007": "Every-third-prime Enots--Wolley",
        "X000008": "Every-fourth-prime Enots--Wolley",
        "X000009": "Square-index-prime-only Enots--Wolley",
        "X000010": "Power-of-two-index-prime-only Enots--Wolley",
        "X000011": "Self-power-index-prime-only Enots--Wolley",
        "X000012": "Full-return Enots--Wolley (2,3)",
        "X000013": "Full-return Enots--Wolley (2,5)",
        "X000014": "Full-return Enots--Wolley (3,5)",
    }

    assert {sequence_id: registry.resolve(sequence_id).name for sequence_id in expected} == expected
