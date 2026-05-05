from abc import ABC, abstractmethod
from dataclasses import dataclass

from .internals.sra import SequenceRead


SRALines = list[tuple[str, ...]]


@dataclass(slots=True)
class SRAOutputUnit(ABC):
    row: SequenceRead

    @property
    def lines(self) -> list[str]:
        return [self.line(i) for i in range(len(self))]

    def __str__(self) -> str:
        return '\n'.join(self.lines) + '\n'

    @abstractmethod
    def __len__(self) -> int:
        ...

    @abstractmethod
    def line(self, idx: int) -> str:
        ...

    @property
    def idx(self) -> int:
        return self.row.idx

    @property
    def split_num(self) -> int:
        return self.row.split_num

    @property
    def header(self) -> str:
        return f'{self.row.accession}.{self.idx} {self.row.name} length={len(self.row.read)}'


class FastaLineTuple(SRAOutputUnit):
    def __len__(self) -> int:
        return 2

    def line(self, idx: int) -> str:
        match idx:
            case 0:
                return f'>{self.header}'
            case 1:
                return self.row.read
            case _:
                msg = "Tuple line index out of range."
                raise IndexError(msg)


class FastqLineQuadruple(SRAOutputUnit):
    def __len__(self) -> int:
        return 4

    def line(self, idx: int) -> str:
        match idx:
            case 0:
                return f'@{self.header}'
            case 1:
                return self.row.read
            case 2:
                return f'+{self.header}'
            case 3:
                return self.row.qual
            case _:
                msg = "Quadruple line index out of range."
                raise IndexError(msg)
