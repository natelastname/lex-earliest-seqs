import pyarrow.parquet as pq

from lex_earliest_seqs import open_run, registry
from lex_earliest_seqs.term_export import write_terms_csv, write_terms_parquet


def test_parquet_export_writes_multiple_batches(tmp_path):
    run = open_run(registry.resolve("ew"), use_cache=False)
    output = tmp_path / "ew.parquet"

    write_terms_parquet(run, 0, 5, output, batch_size=2)

    assert pq.read_table(output).to_pydict() == {
        "subscript": [1, 2, 3, 4, 5],
        "value": [1, 2, 6, 15, 35],
    }


def test_csv_export_respects_sequence_slice(tmp_path):
    run = open_run(registry.resolve("ew"), use_cache=False)
    output = tmp_path / "ew.csv"

    write_terms_csv(run, 2, 5, output, batch_size=2)

    assert output.read_text(encoding="utf-8").splitlines() == [
        "subscript,value",
        "3,6",
        "4,15",
        "5,35",
    ]
