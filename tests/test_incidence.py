import json

from lex_earliest_seqs.core import TermRecord
from lex_earliest_seqs.incidence import (
    ColumnMode,
    build_incidence_table,
    render_json,
    render_text,
)
from lex_earliest_seqs.projections import (
    binary_digit_projection,
    prime_exponent_projection,
)


def records(values):
    return tuple(
        TermRecord(position=i, subscript=i + 1, object_rank=value - 1, value=value)
        for i, value in enumerate(values)
    )


def test_prime_incidence_chronology_records_exponents():
    table = build_incidence_table(
        records([12, 25, 14]),
        projection=prime_exponent_projection(),
    )
    assert table.features == (2, 3, 5, 7)
    assert table.matrix() == (
        (2, 1, None, None),
        (None, None, 2, None),
        (1, None, None, 1),
    )
    rendered = render_text(table)
    assert "object" in rendered
    assert "12" in rendered
    assert "25" in rendered


def test_binary_projection_supports_through_largest_columns():
    table = build_incidence_table(
        records([1, 8]),
        projection=binary_digit_projection(),
        column_mode=ColumnMode.THROUGH_LARGEST,
    )
    assert table.features == (0, 1, 2, 3)


def test_json_preserves_object_rank_metadata():
    table = build_incidence_table(
        records([6]),
        projection=prime_exponent_projection(),
    )
    payload = json.loads(render_json(table))
    assert payload["rows"][0]["object_rank"] == 5
    assert payload["rows"][0]["coordinates"] == {"2": "1", "3": "1"}
