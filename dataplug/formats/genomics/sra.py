from __future__ import annotations

import io
import logging
import os
import ctypes
import tempfile
import vdb

from typing import TYPE_CHECKING
from pathlib import Path
from contextlib import contextmanager
from dataclasses import dataclass
from ctypes import c_ssize_t, c_int, c_void_p, c_size_t, c_uint64, CFUNCTYPE

from ...entities import CloudDataFormat, CloudObjectSlice, PartitioningStrategy
from ...preprocessing.metadata import PreprocessingMetadata

if TYPE_CHECKING:
    from typing import List, Tuple
    from ...cloudobject import CloudObject
#     from botocore.response import StreamingBody

logger = logging.getLogger(__name__)


libvdb_path = os.environ.get('NCBI_VDB_SO_PATH')
assert libvdb_path, (
    'Provide path to ncbi-vdb shared '
    'object in NCBI_VDB_SO_PATH environment variable'
)

manager = vdb.manager(mode=vdb.OpenMode.Read, path=libvdb_path)

default_acc = 'SRR'


class CPPFunctions:
    lib = ctypes.cdll.LoadLibrary(libvdb_path)
    pread_t = CFUNCTYPE(
        c_ssize_t,
        c_int,
        c_void_p,
        c_size_t,
        c_size_t
    )
    mmap_t = CFUNCTYPE(
        c_void_p,
        c_void_p,
        c_size_t,
        c_int,
        c_int,
        c_int,
        c_uint64
    )
    size_t = CFUNCTYPE(c_uint64)

    c_pread = None
    c_size = None
    c_mmap = None

    def __init__(self):
        self.content = None
        self.current = []
        self.analysis = []

        self._register_types()

    def aflush(self):
        if self.current:
            self.analysis.append(self.current)
        self.current = []

    def apread(self, fd, buf, count, offset):
        self.current.append([offset, count])
        return self.pread(fd, buf, count, offset)

    def pread(self, fd, buf, count, offset):
        if self.content is None:
            return -1
        ctypes.memmove(buf, self.content[offset:], count)
        return count

    def mmap(self, addr, length, prot, flags, fd, offset):
        return None

    def size(self):
        return len(self.content)

    def _register_types(self):
        self.lib.register_s3_pread.restype = None
        self.lib.register_s3_pread.argtypes = [self.pread_t]
        self.lib.register_s3_mmap.restype = None
        self.lib.register_s3_mmap.argtypes = [self.mmap_t]

    def _register(self):
        self.lib.register_s3_pread(
            self.c_pread,
            self.c_size
        )
        self.lib.register_s3_mmap(
            self.c_mmap
        )

    def aregister(self):
        self.c_pread = self.pread_t(self.apread)
        self.c_size = self.size_t(self.size)
        self.c_mmap = self.mmap_t(self.mmap)

        self._register()

    def register(self):
        self.c_pread = self.pread_t(self.pread)
        self.c_size = self.size_t(self.size)
        self.c_mmap = self.mmap_t(self.mmap)

        self._register()


@dataclass
class QuadrupleCursor:
    read: vdb.VColumn
    qual: vdb.VColumn
    name: vdb.VColumn

    _column_names = [
        'READ',
        '(INSDC:quality:text:phred_33)QUALITY',
        'NAME',
    ]

    @classmethod
    def from_filepath(cls, filepath):
        table = cls.get_table(filepath)

        columns = table.CreateCursor().OpenColumns(cls._column_names)

        read, qual, name = (columns[name] for name in cls._column_names)

        return cls(read, qual, name)

    @staticmethod
    def get_table(filepath):
        try:
            with manager.OpenDB(filepath) as db:
                return db.OpenTable('SEQUENCE')
        except vdb.vdb_error as e:
            logger.debug(e)
            pass
        try:
            return manager.OpenTable(filepath)
        except vdb.vdb_error as e:
            logger.debug(e)
            pass
        raise IOError(f'Could not open {filepath}')

    def __len__(self):
        return self.read.row_range()[1]

    def row(self, idx, acc=default_acc):
        read = self.read.Read(idx)
        qual = self.qual.Read(idx)
        name = self.name.Read(idx)
        return Quadruple(read=read, qual=qual, name=name, acc=acc)


class QCRange:
    def __init__(self, cursor, start=None, stop=None, step=1, acc=default_acc):
        if start is None:
            start = len(cursor)
        if stop is None:
            start, stop = 0, start

        start += 1
        stop += 1
        self.cursor = cursor
        self.start = start
        self.stop = stop
        self.step = step
        self.current = start
        self.acc = acc

    def __iter__(self):
        return self

    def _finished(self):
        return (
            (self.step > 0 and self.current >= self.stop)
            or
            (self.step < 0 and self.current <= self.stop)
        )

    def __next__(self):
        if self._finished():
            raise StopIteration

        result = self.cursor.row(self.current, acc=self.acc)
        self.current += self.step
        return result


@dataclass
class Quadruple:
    read: str
    qual: str
    name: str
    acc: str = default_acc

    def __str__(self):
        return '\n'.join(self.lines)

#     def __str__(self):
#         return f'''@{self.acc}.{self.name} length={len(self.read)}
# {self.read}
# +{self.acc}.{self.name} length={len(self.read)}
# {self.qual}'''

    @property
    def lines(self) -> List[str]:
        return [
            self.line0,
            self.line1,
            self.line2,
            self.line3,
        ]

    @property
    def line0(self):
        return f'@{self.acc}.{self.name} length={len(self.read)}'

    @property
    def line1(self):
        return f'{self.read}'

    @property
    def line2(self):
        return f'+{self.acc}.{self.name} length={len(self.read)}'

    @property
    def line3(self):
        return f'{self.qual}'


def fasterq_dump(filepath):
    c = QuadrupleCursor.from_filepath(filepath)

    for q in QCRange(c, 12):
        pass
        print(q)


@contextmanager
def temporary_sra():
    with tempfile.TemporaryDirectory() as dpath:
        fpath = Path(dpath) / 'SRR21370000'
        fp = open(fpath, 'w+')
        fp.write('s3re')
        fp.flush()
        yield str(fpath)


def download(cloud_object: CloudObject) -> bytes:
    bucket = cloud_object.path.bucket
    key = cloud_object.path.key
    obj_resp = cloud_object.storage.get_object(Bucket=bucket, Key=key)
    assert obj_resp.get("ResponseMetadata", {}).get("HTTPStatusCode") == 200
    data_stream = obj_resp["Body"]

    CHUNK_SIZE = 655360000
    content = b""
    chunk = data_stream.read(CHUNK_SIZE)
    while chunk != b"":
        content += chunk
        chunk = data_stream.read(CHUNK_SIZE)
        logger.info(f'Next bytes {len(content)}')
    if hasattr(data_stream, "close"):
        data_stream.close()

    return content


def preprocess_sra(cloud_object: CloudObject):
    logger.info('Preprocessing sra started')

    functions = CPPFunctions()
    functions.content = download(cloud_object)
    functions.aregister()

    import tqdm
    with temporary_sra() as fpath:
        qc = QuadrupleCursor.from_filepath(fpath)
        total_lines = len(qc)
        functions.aflush()
        for row in tqdm.tqdm(QCRange(qc, total_lines), total=total_lines):
            functions.aflush()
    import pprint
    pprint.pprint(functions.analysis)
    print(f'Analysis length: {len(functions.analysis)}')
    return PreprocessingMetadata(
        metadata=io.BytesIO(b'nempty'),
        attributes={
            'total_lines': total_lines,
        }
    )


@CloudDataFormat(preprocessing_function=preprocess_sra)
class SRA:
    total_lines: int
    index_key: str


class SRASlice(CloudObjectSlice):
    def __init__(self, start, end, *args, **kwargs):
        self.start = start
        self.end = end
        super().__init__(*args, **kwargs)

    def get(self) -> list[str]:
        functions = CPPFunctions()
        functions.content = download(self.cloud_object)
        functions.register()

        lines = []
        with temporary_sra() as fpath:
            qc = QuadrupleCursor.from_filepath(fpath)
            for row in QCRange(qc, self.start, self.end):
                lines.extend(row.lines)

        return lines


def partition_into_ranges(
    total_lines: int,
    num_chunks: int
) -> List[Tuple[int, int]]:
    n, r = divmod(total_lines, num_chunks)

    if total_lines == 0:
        return []
    if total_lines < num_chunks:
        logger.warning(
            'Number of chunks given is greater '
            'than total number of lines in file. '
            'Dividing files into single-line chunks... '
            'For performance gains, consider reducing '
            'number of chunks.'
        )
        num_chunks = total_lines

    n, r = divmod(total_lines, num_chunks)  # lines_per_chunk, leftovers
    firsts = [
        ((n + 1) * i, (n + 1) * (i + 1))
        for i in range(r)
    ]
    last = firsts[-1][-1] if firsts else 0
    seconds = [
        (last + n * i, last + n * (i + 1))
        for i in range(num_chunks - r)
    ]
    return firsts + seconds


@PartitioningStrategy(dataformat=SRA)
def partition_chunks_strategy(
    cloud_object: CloudObject,
    num_chunks: int
) -> List[SRASlice]:
    logger.info('SRA partitioning started')
    total_lines = int(cloud_object.get_attribute("total_lines"))
    ranges = partition_into_ranges(total_lines, num_chunks)
    return [SRASlice(*range) for range in ranges]
