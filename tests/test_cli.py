import pyarrow.parquet as pq
import pytest

from lex_earliest_seqs.cli import app


def _run_cli(args: list[str]) -> None:
    """Invoke Cyclopts as a CLI and require a successful process exit."""

    with pytest.raises(SystemExit) as exc_info:
        app(args)
    assert exc_info.value.code == 0


def test_list_command(capsys):
    _run_cli(["list"])

    output = capsys.readouterr().out
    assert "A336957" in output
    assert "A338833" in output


def test_terms_command_without_cache(capsys):
    _run_cli(["terms", "ew", "5", "--no-cache"])

    captured = capsys.readouterr()
    assert captured.out.splitlines() == [
        "1\t1",
        "2\t2",
        "3\t6",
        "4\t15",
        "5\t35",
    ]
    assert "compute A336957:" in captured.err
    assert "5/5 terms (100.0%)" in captured.err


def test_progress_can_be_disabled(capsys):
    _run_cli(["terms", "ew", "5", "--no-cache", "--no-progress"])

    captured = capsys.readouterr()
    assert captured.err == ""


def test_terms_csv_output_with_short_output_option(tmp_path, capsys):
    output = tmp_path / "ew.csv"
    _run_cli(
        [
            "terms",
            "ew",
            "5",
            "--no-cache",
            "--no-progress",
            "-o",
            str(output),
            "--format",
            "csv",
        ]
    )

    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""
    assert output.read_text(encoding="utf-8").splitlines() == [
        "subscript,value",
        "1,1",
        "2,2",
        "3,6",
        "4,15",
        "5,35",
    ]


def test_terms_output_format_is_inferred_from_csv_suffix(tmp_path):
    output = tmp_path / "ew.csv"
    _run_cli(
        [
            "terms",
            "ew",
            "3",
            "--no-cache",
            "--no-progress",
            "--output",
            str(output),
        ]
    )

    assert output.read_text(encoding="utf-8").splitlines() == [
        "subscript,value",
        "1,1",
        "2,2",
        "3,6",
    ]


def test_terms_parquet_output_preserves_subscript_and_value(tmp_path, capsys):
    output = tmp_path / "ew.data"
    _run_cli(
        [
            "terms",
            "ew",
            "3",
            "--start-position",
            "2",
            "--no-cache",
            "--no-progress",
            "-o",
            str(output),
            "--format",
            "parquet",
        ]
    )

    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""
    table = pq.read_table(output)
    assert table.column_names == ["subscript", "value"]
    assert table.to_pydict() == {
        "subscript": [3, 4, 5],
        "value": [6, 15, 35],
    }


def test_terms_output_format_is_inferred_from_parquet_suffix(tmp_path):
    output = tmp_path / "ew.parquet"
    _run_cli(
        [
            "terms",
            "ew",
            "2",
            "--no-cache",
            "--no-progress",
            "-o",
            str(output),
        ]
    )

    assert pq.read_table(output).to_pydict() == {
        "subscript": [1, 2],
        "value": [1, 2],
    }
