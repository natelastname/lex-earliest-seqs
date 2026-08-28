from dataclasses import dataclass, field

import pytest

from lex_earliest_seqs.core import SequenceDefinition, SequenceRegistry, SequenceRun
from lex_earliest_seqs.object_space import PositiveIntegers


@dataclass
class CountingGenerator:
    terms: list[int] = field(default_factory=list)

    def extend_to(self, count: int) -> None:
        while len(self.terms) < count:
            self.terms.append(len(self.terms) + 1)


def definition() -> SequenceDefinition[int]:
    return SequenceDefinition(
        id="counting",
        name="Counting",
        aliases=("count",),
        generator_factory=CountingGenerator,
        object_space=PositiveIntegers(),
        offset=1,
    )


def test_position_subscript_and_object_rank_are_distinct_fields():
    run = SequenceRun(definition(), CountingGenerator())
    record = run.record_at_position(4)
    assert record.position == 4
    assert record.subscript == 5
    assert record.object_rank == 4
    assert record.value == 5


def test_registry_resolves_id_name_and_alias():
    sequence_registry = SequenceRegistry()
    sequence_registry.register(definition())
    assert sequence_registry.resolve("COUNTING").id == "counting"
    assert sequence_registry.resolve("Counting").id == "counting"
    assert sequence_registry.resolve("count").id == "counting"


def test_invalid_subscript_fails():
    run = SequenceRun(definition(), CountingGenerator())
    with pytest.raises(IndexError):
        run.at_subscript(0)
