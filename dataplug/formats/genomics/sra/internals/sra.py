from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Any, Protocol, TypeGuard, TypeVar, override

from ._loader import get_table

VColumnRead = Callable[[int], Any]
T = TypeVar("T")


def no_nones(items: Iterable[T | None]) -> TypeGuard[Iterable[T]]:
    return all(item is not None for item in items)


@dataclass(slots=True)
class SequenceRead:
    _spot: Spot
    start: int
    length: int
    split_num: int

    @property
    def read(self) -> str:
        return self._spot.read[self.start:self.start + self.length]

    @property
    def qual(self) -> str:
        return self._spot.qual[self.start:self.start + self.length]

    @property
    def name(self) -> str:
        return self._spot.name

    @property
    def accession(self) -> str:
        return self._spot.accession

    @property
    def idx(self) -> int:
        return self._spot.idx


class VColumn(Protocol):
    def Read(self, idx: int) -> Any:  # noqa: N802
        ...

    def row_range(self) -> tuple[int, int]:
        ...


class FakeVColumn(VColumn):
    def __init__(self, read_func: VColumnRead):
        self._read_func = read_func

    def Read(self, idx: int) -> Any:  # noqa: N802
        return self._read_func(idx)

    @override
    def row_range(self) -> tuple[int, int]:
        return (1, sys.maxsize)


@dataclass(slots=True)
class Spot:
    accession: str
    idx: int
    read: str
    qual: str
    name: str
    starts: list[int]
    lengths: list[int]

    def reads(self) -> list[SequenceRead]:
        reads = []
        for i, (start, length) in enumerate(zip(self.starts, self.lengths, strict=False)):
            if length > 0:
                reads.append(SequenceRead(
                    _spot=self,
                    start=start,
                    length=length,
                    split_num=i + 1
                ))

        return reads


class SequenceTableVCursor:
    column_names = (
        "READ",
        "(INSDC:quality:text:phred_33)QUALITY",
        "NAME",
        "READ_START",
        "READ_LEN"
        )

    def __init__(self, accession: str, split: bool = True):
        self.accession = accession
        self.split = split
        self.columns = self._unpack_columns()

    def _columns(self) -> dict[str, VColumn]:
        table = get_table(self.accession)
        cursor = table.CreateCursor()
        return cursor.OpenColumns(list(self.column_names))

    def _unpack_columns(self) -> tuple[VColumn, ...]:
        columns_dict = self._columns()
        columns = [columns_dict.get(name) for name in self.column_names]

        if not self.split or columns[-1] is None or columns[-2] is None:
            columns[-2] = FakeVColumn(self._default_read_starts)
            columns[-1] = FakeVColumn(self._default_read_lengths)

        if not no_nones(columns):
            msg = (
                "Could not open sequence table. "
                f"{self.accession} is not an SRA-object"
            )
            raise ValueError(msg)

        return tuple(columns)

    def spot(self, idx: int) -> Spot:
        data = [col.Read(idx) for col in self.columns]
        return Spot(self.accession, idx, *data)

    @staticmethod
    def _default_read_starts(_: int) -> list[int]:
        return [0]

    def _default_read_lengths(self, idx: int) -> list[int]:
        return [len(self.read_column.Read(idx))]

    @property
    def read_column(self) -> VColumn:
        return self.columns[0]

    def __len__(self):
        return self.read_column.row_range()[1]


def _init_read_range(*args: int) -> range:
    """ Initialize a range for read indices based on arguments.
        Mimics Python's built-in range.
    """
    argc = len(args)

    if argc == 0:
        msg = 'QuadrupleRange expected at least 1 argument, got 0'
        raise TypeError(msg)
    if argc > 3:  # noqa: PLR2004
        msg = f'QuadrupleRange expected at most 3 arguments, got {len(args)}'
        raise TypeError(msg)

    for i, arg in enumerate(args):
        if not isinstance(arg, int):
            msg = f'QuadrupleRange argument {i} must be an integer'
            raise TypeError(msg)

    match argc:
        case 1:
            start = 1
            stop = args[0] + 1
            step = 1
        case 2:
            start = args[0] + 1
            stop = args[1] + 1
            step = 1
        case 3:
            start = args[0] + 1
            stop = args[1] + 1
            step = args[2]
        case _:
            msg = "Unreachable code"
            raise RuntimeError(msg)

    if step == 0:
        msg = "QuadrupleRange argument 3 must not be zero"
        raise ValueError()
    return range(start, stop, step)


class SequenceTableVRange:
    def __init__(self, cursor: SequenceTableVCursor, *args: int):
        self.cursor = cursor
        self.range = _init_read_range(*args)

    def __len__(self):
        return len(self.range)

    def __iter__(self):
        for idx in self.range:
            yield self.cursor.spot(idx)


class SRAFile:
    def __init__(self, filepath: str | Path, paired: bool = True):
        accession = Path(filepath).name
        self.accession = accession
        self.paired = paired
        self._cursor: SequenceTableVCursor | None = None

    def range(self, *args: int) -> SequenceTableVRange:
        if self._cursor is None:
            self._cursor = SequenceTableVCursor(self.accession, split=self.paired)
        return SequenceTableVRange(self._cursor, *args)

    def __len__(self):
        if self._cursor is None:
            self._cursor = SequenceTableVCursor(self.accession, split=self.paired)
        return len(self._cursor)
