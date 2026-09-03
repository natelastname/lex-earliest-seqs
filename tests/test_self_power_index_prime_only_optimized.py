from lex_earliest_seqs import registry
from lex_earliest_seqs.zoo import sparse_prime_index_only_enots_wolley as sparse_module
from lex_earliest_seqs.zoo.self_power_index_prime_only_optimized import (
    SELF_POWER_RETAINED_PRIME_PREFIX,
    SelfPowerIndexPrimeOnlyEnotsWolleyGenerator,
)


def test_precomputed_self_power_prime_coordinates_are_used_without_nth_prime(
    monkeypatch,
):
    def fail_nth_prime(_index):
        raise AssertionError("precomputed self-power coordinate called nth_prime")

    monkeypatch.setattr(sparse_module, "nth_prime", fail_nth_prime)
    generator = SelfPowerIndexPrimeOnlyEnotsWolleyGenerator()

    actual = [
        generator._prime_at_position(position)
        for position in range(len(SELF_POWER_RETAINED_PRIME_PREFIX))
    ]

    assert actual == list(SELF_POWER_RETAINED_PRIME_PREFIX)
    assert [generator._allowed_prime_index(i) for i in range(10)] == [
        1,
        4,
        27,
        256,
        3_125,
        46_656,
        823_543,
        16_777_216,
        387_420_489,
        10_000_000_000,
    ]


def test_exact_precomputed_prime_provides_candidate_cutoff_without_materializing():
    generator = SelfPowerIndexPrimeOnlyEnotsWolleyGenerator()

    assert generator._prime_at_position(8, upper_bound=8_448_283_756) is None
    assert generator.retained_primes == []

    assert generator._prime_at_position(8, upper_bound=8_448_283_757) == 8_448_283_757
    assert generator.retained_primes == list(SELF_POWER_RETAINED_PRIME_PREFIX[:9])


def test_registered_x000014_uses_cached_self_power_generator():
    definition = registry.resolve("X000014")
    generator = definition.generator_factory()

    assert type(generator) is SelfPowerIndexPrimeOnlyEnotsWolleyGenerator
    assert definition.generator_version == 3


def test_cached_self_power_generator_preserves_known_prefix():
    generator = SelfPowerIndexPrimeOnlyEnotsWolleyGenerator()
    generator.extend_to(20)

    assert generator.terms == [
        1,
        2,
        14,
        721,
        166_757,
        3_238,
        28,
        5_047,
        2_954_761,
        57_374,
        56,
        11_333,
        17_175_971,
        206,
        98,
        79_331,
        46_444_253,
        114_748,
        112,
        35_329,
    ]
