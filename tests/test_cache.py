from dataclasses import dataclass, field

import pytest

from lex_earliest_seqs.cache import CacheCompatibilityError, open_run
from lex_earliest_seqs.core import SequenceDefinition
from lex_earliest_seqs.object_space import PositiveIntegers


@dataclass
class HugeIntegerGenerator:
    terms: list[int] = field(default_factory=list)
    next_exponent: int = 100

    def extend_to(self, count: int) -> None:
        while len(self.terms) < count:
            self.terms.append(2**self.next_exponent)
            self.next_exponent += 100


def definition(*, generator_version: int = 1) -> SequenceDefinition[int]:
    return SequenceDefinition(
        id="huge-integers",
        name="Huge integers",
        generator_factory=HugeIntegerGenerator,
        generator_version=generator_version,
        object_space=PositiveIntegers(),
    )


def test_pickle_cache_preserves_generator_state_and_arbitrary_precision(tmp_path):
    first = open_run(definition(), cache_dir=tmp_path)
    first.ensure(3)
    assert first.terms[-1] == 2**300

    second = open_run(definition(), cache_dir=tmp_path)
    assert list(second.terms) == [2**100, 2**200, 2**300]
    second.ensure(4)
    assert second.terms[-1] == 2**400
    assert second.generator.next_exponent == 500


def test_generator_version_mismatch_fails_loudly(tmp_path):
    run = open_run(definition(generator_version=1), cache_dir=tmp_path)
    run.ensure(1)
    with pytest.raises(CacheCompatibilityError):
        open_run(definition(generator_version=2), cache_dir=tmp_path)


def test_generated_prefix_is_written(tmp_path):
    run = open_run(definition(), cache_dir=tmp_path)
    run.ensure(1)
    assert run.cache_path is not None and run.cache_path.exists()


def test_pickle_load_reports_byte_progress(tmp_path):
    first = open_run(definition(), cache_dir=tmp_path)
    first.ensure(4)

    events: list[tuple[int, int]] = []
    second = open_run(
        definition(),
        cache_dir=tmp_path,
        load_progress=lambda current, total: events.append((current, total)),
    )

    assert second.terms[-1] == 2**400
    assert events[0][0] == 0
    assert events[-1][0] == events[-1][1]
    assert events[-1][1] > 0
    assert all(0 <= current <= total for current, total in events)
