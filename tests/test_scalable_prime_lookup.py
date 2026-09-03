from lex_earliest_seqs.zoo import scalable_prime_lookup


def test_segmented_fallback_matches_known_nth_primes(monkeypatch):
    monkeypatch.setenv("LEX_EARLIEST_SEQS_DISABLE_PRIMECOUNT", "1")
    monkeypatch.setenv("LEX_EARLIEST_SEQS_DISABLE_PRIMESIEVE", "1")
    scalable_prime_lookup._reset_for_tests()

    expected = {
        1: 2,
        10: 29,
        100: 541,
        1_000: 7_919,
        10_000: 104_729,
        100_000: 1_299_709,
    }
    for index, prime in expected.items():
        assert scalable_prime_lookup.scalable_nth_prime(index) == prime


def test_segmented_fallback_handles_out_of_order_queries(monkeypatch):
    monkeypatch.setenv("LEX_EARLIEST_SEQS_DISABLE_PRIMECOUNT", "1")
    monkeypatch.setenv("LEX_EARLIEST_SEQS_DISABLE_PRIMESIEVE", "1")
    scalable_prime_lookup._reset_for_tests()

    assert scalable_prime_lookup.scalable_nth_prime(10_000) == 104_729
    assert scalable_prime_lookup.scalable_nth_prime(1_000) == 7_919
    assert scalable_prime_lookup.scalable_nth_prime(20_000) == 224_737


def test_optional_primesieve_backend_advances_from_nearest_known_prime(monkeypatch):
    calls: list[tuple[int, int]] = []

    class FakePrimeSieve:
        @staticmethod
        def nth_prime(index_gap: int, start: int = 0) -> int:
            calls.append((index_gap, start))
            values = {
                (9, 2): 29,
                (90, 29): 541,
                (900, 541): 7_919,
            }
            return values[(index_gap, start)]

    scalable_prime_lookup._reset_for_tests()
    monkeypatch.setenv("LEX_EARLIEST_SEQS_DISABLE_PRIMECOUNT", "1")
    monkeypatch.delenv("LEX_EARLIEST_SEQS_DISABLE_PRIMESIEVE", raising=False)
    monkeypatch.setattr(
        scalable_prime_lookup,
        "_load_primesieve",
        lambda: FakePrimeSieve,
    )

    assert scalable_prime_lookup.scalable_nth_prime(10) == 29
    assert scalable_prime_lookup.scalable_nth_prime(100) == 541
    assert scalable_prime_lookup.scalable_nth_prime(1_000) == 7_919
    assert calls == [(9, 2), (90, 29), (900, 541)]


def test_primecount_backend_wins_for_huge_indices(monkeypatch):
    scalable_prime_lookup._reset_for_tests()
    monkeypatch.delenv("LEX_EARLIEST_SEQS_DISABLE_PRIMECOUNT", raising=False)
    monkeypatch.setattr(
        scalable_prime_lookup,
        "_load_primecount_executable",
        lambda: "/fake/primecount",
    )

    calls: list[tuple[str, int]] = []

    def fake_primecount(executable: str, index: int) -> int:
        calls.append((executable, index))
        return 15_485_863

    monkeypatch.setattr(
        scalable_prime_lookup,
        "_primecount_nth_prime",
        fake_primecount,
    )
    monkeypatch.setattr(
        scalable_prime_lookup,
        "_load_primesieve",
        lambda: (_ for _ in ()).throw(
            AssertionError("primesieve should not run for huge primecount query")
        ),
    )

    index = scalable_prime_lookup._PRIMECOUNT_INDEX_THRESHOLD
    assert scalable_prime_lookup.scalable_nth_prime(index) == 15_485_863
    assert calls == [("/fake/primecount", index)]


def test_primecount_output_parser(monkeypatch):
    class Completed:
        stdout = "15,485,863\n"

    monkeypatch.setattr(scalable_prime_lookup, "run", lambda *args, **kwargs: Completed())
    assert scalable_prime_lookup._primecount_nth_prime("primecount", 1_000_000) == 15_485_863
