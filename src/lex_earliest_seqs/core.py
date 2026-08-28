"""Core abstractions for lexicographically earliest sequences."""

from __future__ import annotations

from collections.abc import Callable, Iterator, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Generic, Protocol, TypeVar, runtime_checkable

ObjectT = TypeVar("ObjectT")


@runtime_checkable
class SequenceGenerator(Protocol[ObjectT]):
    """Mutable, pickleable state for computing a sequence prefix."""

    terms: list[ObjectT]

    def extend_to(self, count: int) -> None:
        """Extend ``terms`` so it contains at least ``count`` terms."""


@runtime_checkable
class ObjectSpace(Protocol[ObjectT]):
    """Ordered ambient space of admissible objects."""

    key: str

    def at_rank(self, rank: int) -> ObjectT:
        """Return the object having zero-based ambient ``rank``."""

    def rank_of(self, value: ObjectT) -> int:
        """Return the zero-based ambient rank of ``value``."""


@dataclass(frozen=True, slots=True)
class TermRecord(Generic[ObjectT]):
    """A sequence occurrence with both sequence and ambient coordinates."""

    position: int
    subscript: int
    object_rank: int
    value: ObjectT


@dataclass(frozen=True)
class SequenceDefinition(Generic[ObjectT]):
    """Stable mathematical metadata plus factories and display projections."""

    id: str
    name: str
    generator_factory: Callable[[], SequenceGenerator[ObjectT]]
    object_space: ObjectSpace[ObjectT]
    generator_version: int = 1
    definition_version: int = 1
    offset: int = 1
    oeis: str | None = None
    aliases: tuple[str, ...] = ()
    projections: Mapping[str, object] | None = None
    description: str = ""

    def __post_init__(self) -> None:
        if not self.id.strip():
            raise ValueError("sequence id must be non-empty")
        if not self.name.strip():
            raise ValueError("sequence name must be non-empty")
        if self.generator_version < 1:
            raise ValueError("generator_version must be positive")
        if self.definition_version < 1:
            raise ValueError("definition_version must be positive")
        object.__setattr__(self, "aliases", tuple(self.aliases))
        object.__setattr__(self, "projections", dict(self.projections or {}))


class SequenceRegistry:
    """Registry of immutable sequence definitions and their aliases."""

    def __init__(self) -> None:
        self._definitions: dict[str, SequenceDefinition[object]] = {}
        self._aliases: dict[str, str] = {}

    @staticmethod
    def _normalize(name: str) -> str:
        return name.strip().lower()

    def register(self, definition: SequenceDefinition[object]) -> None:
        canonical = self._normalize(definition.id)
        names = [definition.id, definition.name, *definition.aliases]
        if definition.oeis:
            names.append(definition.oeis)
        normalized = [self._normalize(name) for name in names]

        if canonical in self._definitions:
            raise ValueError(f"sequence already registered: {definition.id}")
        collisions = [name for name in normalized if name in self._aliases]
        if collisions:
            raise ValueError(f"sequence alias already registered: {collisions[0]}")

        self._definitions[canonical] = definition
        for name in normalized:
            self._aliases[name] = canonical

    def resolve(self, name: str) -> SequenceDefinition[object]:
        normalized = self._normalize(name)
        try:
            canonical = self._aliases[normalized]
        except KeyError as exc:
            raise KeyError(f"unknown sequence: {name}") from exc
        return self._definitions[canonical]

    def definitions(self) -> tuple[SequenceDefinition[object], ...]:
        return tuple(sorted(self._definitions.values(), key=lambda item: item.id))

    def __iter__(self) -> Iterator[SequenceDefinition[object]]:
        return iter(self.definitions())


class SequenceRun(Generic[ObjectT]):
    """A loaded generator plus its definition and optional persistent cache."""

    def __init__(
        self,
        definition: SequenceDefinition[ObjectT],
        generator: SequenceGenerator[ObjectT],
        *,
        cache_path: Path | None = None,
    ) -> None:
        self.definition = definition
        self.generator = generator
        self.cache_path = cache_path

    @property
    def terms(self) -> Sequence[ObjectT]:
        return self.generator.terms

    def ensure(self, count: int, *, save: bool = True) -> Sequence[ObjectT]:
        if count < 0:
            raise ValueError("count must be nonnegative")
        before = len(self.generator.terms)
        self.generator.extend_to(count)
        if len(self.generator.terms) < count:
            raise RuntimeError(
                f"generator produced only {len(self.generator.terms)} terms "
                f"after extend_to({count})"
            )
        if save and self.cache_path is not None and len(self.generator.terms) != before:
            from .cache import save_generator

            save_generator(self.definition, self.generator, self.cache_path)
        return self.generator.terms[:count]

    def at_position(self, position: int) -> ObjectT:
        if position < 0:
            raise IndexError("position must be nonnegative")
        self.ensure(position + 1)
        return self.generator.terms[position]

    def at_subscript(self, subscript: int) -> ObjectT:
        position = subscript - self.definition.offset
        if position < 0:
            raise IndexError(
                f"subscript {subscript} precedes sequence offset {self.definition.offset}"
            )
        return self.at_position(position)

    def record_at_position(self, position: int) -> TermRecord[ObjectT]:
        value = self.at_position(position)
        return TermRecord(
            position=position,
            subscript=position + self.definition.offset,
            object_rank=self.definition.object_space.rank_of(value),
            value=value,
        )

    def records(
        self, start: int = 0, stop: int | None = None
    ) -> tuple[TermRecord[ObjectT], ...]:
        if start < 0:
            raise ValueError("start must be nonnegative")
        if stop is None:
            stop = len(self.generator.terms)
        if stop < start:
            raise ValueError("stop must be at least start")
        self.ensure(stop)
        return tuple(self.record_at_position(position) for position in range(start, stop))
