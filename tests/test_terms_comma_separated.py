import pytest

from lex_earliest_seqs.cli import app


def _run_cli(args: list[str]) -> None:
    """Invoke Cyclopts as a CLI and require a successful process exit."""

    with pytest.raises(SystemExit) as exc_info:
        app(args)
    assert exc_info.value.code == 0


def test_terms_comma_separated_outputs_values_on_one_line(capsys):
    _run_cli(
        [
            "terms",
            "ew",
            "3",
            "--start-position",
            "2",
            "--comma-separated",
            "--no-cache",
            "--no-progress",
        ]
    )

    captured = capsys.readouterr()
    assert captured.out == "6,15,35\n"
    assert captured.err == ""
