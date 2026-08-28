from lex_earliest_seqs.cli import app


def test_list_command(capsys):
    app(["list"])

    output = capsys.readouterr().out
    assert "A336957" in output
    assert "A338833" in output


def test_terms_command_without_cache(capsys):
    app(["terms", "ew", "5", "--no-cache"])

    assert capsys.readouterr().out.splitlines() == [
        "1\t1",
        "2\t2",
        "3\t6",
        "4\t15",
        "5\t35",
    ]
