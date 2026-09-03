import pytest

from lex_earliest_seqs.zoo import scalable_prime_lookup


def test_primecountpy_backend_matches_known_nth_primes():
    scalable_prime_lookup._reset_for_tests()

    expected = {
        1: 2,
        10: 29,
        100: 541,
        1_000: 7_919,
        100_000: 1_299_709,
        1_000_000: 15_485_863,
        387_420_489: 8_448_283_757,
    }
    for index, prime in expected.items():
        assert scalable_prime_lookup.scalable_nth_prime(index) == prime


def test_scalable_nth_prime_delegates_to_primecountpy_and_caches(monkeypatch):
    calls: list[int] = []

    def fake_nth_prime(index: int) -> int:
        calls.append(index)
        return 15_485_863

    scalable_prime_lookup._reset_for_tests()
    monkeypatch.setattr(scalable_prime_lookup, "primecount_nth_prime", fake_nth_prime)

    assert scalable_prime_lookup.scalable_nth_prime(1_000_000) == 15_485_863
    assert scalable_prime_lookup.scalable_nth_prime(1_000_000) == 15_485_863
    assert calls == [1_000_000]


def test_scalable_nth_prime_rejects_invalid_indices():
    scalable_prime_lookup._reset_for_tests()

    with pytest.raises(TypeError):
        scalable_prime_lookup.scalable_nth_prime(3.0)
    with pytest.raises(ValueError):
        scalable_prime_lookup.scalable_nth_prime(0)


def test_scalable_nth_prime_rejects_invalid_backend_result(monkeypatch):
    scalable_prime_lookup._reset_for_tests()
    monkeypatch.setattr(scalable_prime_lookup, "primecount_nth_prime", lambda _index: -1)

    with pytest.raises(RuntimeError, match="invalid nth prime"):
        scalable_prime_lookup.scalable_nth_prime(1_000_000)
