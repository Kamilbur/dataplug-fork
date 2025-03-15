from __future__ import annotations

import io
import logging
import os
import sys
import ctypes
import tempfile
import struct
import vdb

from typing import TYPE_CHECKING
from pathlib import Path
from contextlib import contextmanager
from dataclasses import dataclass
from ctypes import c_int, c_uint64, c_void_p, c_size_t, c_char_p, Structure, pointer, POINTER, sizeof
from collections import defaultdict

from ...entities import CloudDataFormat, CloudObjectSlice, PartitioningStrategy
from ...preprocessing.metadata import PreprocessingMetadata
from multiprocessing.managers import SharedMemoryManager

if TYPE_CHECKING:
    from typing import List, Tuple
    from ...cloudobject import CloudObject

logger = logging.getLogger(__name__)

CHUNK_SIZE = 1048576
PREAD_BUF_SIZE = 100 * sizeof(c_uint64)
MMAP_BUF_SIZE = 100 * sizeof(c_uint64)

libvdb_path = os.environ.get('NCBI_VDB_SO_PATH')
assert libvdb_path, (
    'Provide path to ncbi-vdb shared '
    'object in NCBI_VDB_SO_PATH environment variable'
)

default_acc = 'SRR'


_endianness = '<' if sys.byteorder == 'little' else '>'
_uint64_t_fmt = _endianness + 'Q'
_size_t_fmt = _endianness
if sys.maxsize < 2**32:
    _size_t_fmt += 'I'
else:
    _size_t_fmt += 'Q'


_ncbi_vdb_mapping = None


class ShmInfo(Structure):
    _fields_ = [
        ("fd", c_int),
        ("ptr", c_void_p),
        ("length", c_size_t),
    ]


ShmInfo_p = POINTER(ShmInfo)


class NcbiVdbMapping:
    def __init__(self):
        self.lib = ctypes.cdll.LoadLibrary(libvdb_path)
        self.lib.register_shmem.restype = None
        self.lib.register_shmem.argtypes = [ShmInfo_p, c_char_p]

        self.lib.close_shmem.restype = None
        self.lib.close_shmem.argtypes = [ShmInfo_p]

        self._dp_sra_size = c_size_t.in_dll(self.lib, "dp_sra_size")
        self._dp_mode = c_int.in_dll(self.lib, "dp_mode")

        self.shm_buf = ShmInfo.in_dll(self.lib, "shm_buf")
        self.mmap_buf = ShmInfo.in_dll(self.lib, "mmap_buf")
        self.pread_buf = ShmInfo.in_dll(self.lib, "pread_buf")

    @property
    def dp_sra_size(self):
        return self._dp_sra_size.value

    @dp_sra_size.setter
    def dp_sra_size(self, value: int):
        self._dp_sra_size.value = value

    @property
    def dp_mode(self):
        return self._dp_mode.value

    @dp_mode.setter
    def dp_mode(self, value: int):
        self._dp_mode.value = value

    def register_shmem(self, shm_p: ShmInfo_p, name: bytes) -> None:
        self.lib.register_shmem(shm_p, name.encode())

    def close_shmem(self, shm_p: ShmInfo_p) -> None:
        self.lib.close_shmem(shm_p)


def get_ncbi_vdb_mapping():
    global _ncbi_vdb_mapping
    if _ncbi_vdb_mapping is None:
        _ncbi_vdb_mapping = NcbiVdbMapping()
    return _ncbi_vdb_mapping


class NcbiVdbSharedMemory:
    def __init__(self, shm_info: ShmInfo, size: int):
        self.size = size
        self.shm_info = shm_info
        self.shared_memory_manager = SharedMemoryManager()
        self.mapping = get_ncbi_vdb_mapping()

    @property
    def buf(self):
        return self.shared_memory.buf

    @property
    def name(self):
        return self.shared_memory.name

    def __enter__(self):
        self._init_shared_memory()
        self._init_shm_info()
        return self.shared_memory

    def _init_shared_memory(self):
        self.shared_memory_manager.start()
        self.shared_memory = self.shared_memory_manager.SharedMemory(size=self.size)

    def _init_shm_info(self):
        self.shm_info.length = self.size
        self.mapping.register_shmem(pointer(self.shm_info), self.name)

    def __exit__(self, type, value, traceback):
        self.mapping.close_shmem(pointer(self.shm_info))
        self.shared_memory_manager.shutdown()


@dataclass
class QuadrupleCursor:
    read: vdb.VColumn
    qual: vdb.VColumn
    name: vdb.VColumn

    manager = None

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

    @classmethod
    def get_table(cls, filepath):
        if cls.manager is None:
            cls.manager = vdb.manager(mode=vdb.OpenMode.Read, path=libvdb_path)
        try:
            with cls.manager.OpenDB(filepath) as db:
                return db.OpenTable('SEQUENCE')
        except vdb.vdb_error as e:
            logger.debug(e)
            pass
        try:
            return cls.manager.OpenTable(filepath)
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


@contextmanager
def temporary_sra():
    with tempfile.TemporaryDirectory() as dpath:
        fpath = Path(dpath) / 'SRR21370000'
        fp = open(fpath, 'w+')
        fp.write('s3re')
        fp.flush()
        yield str(fpath)


def download(cloud_object: CloudObject, destination: memoryview) -> bytes:
    bucket = cloud_object.path.bucket
    key = cloud_object.path.key
    obj_resp = cloud_object.storage.get_object(Bucket=bucket, Key=key)
    assert obj_resp.get("ResponseMetadata", {}).get("HTTPStatusCode") == 200
    data_stream = obj_resp["Body"]

    chunk = data_stream.read(CHUNK_SIZE)
    idx = 0
    while chunk != b"":
        destination[idx:idx+len(chunk)] = chunk
        idx += len(chunk)
        chunk = data_stream.read(CHUNK_SIZE)
        logger.info(f'Next bytes {len(chunk)}')
    if hasattr(data_stream, "close"):
        data_stream.close()


def unpack_from(buf, pos, size=sizeof(c_uint64), fmt=_uint64_t_fmt):
    return struct.unpack(fmt, buf[pos:pos + size])[0]


def pack_to(buf, pos, val, size=sizeof(c_uint64), fmt=_uint64_t_fmt):
    buf[pos:pos + size] = struct.pack(fmt, val)


def ibuf_to_list(buf, length, total_size):
    retval = []

    idx = 0
    for i in range(length // 2):
        pos = unpack_from(buf, idx)
        size = unpack_from(
            buf,
            idx + sizeof(c_uint64),
            size=sizeof(c_size_t),
            fmt=_size_t_fmt
        )
        # Additional 4 bytes for alignment purposes
        if pos + size + 4 < total_size:
            size += 4
        retval.append((pos, pos + size))
        idx += sizeof(c_uint64) + sizeof(c_size_t)

    return retval


def sum_intervals(intervals):
    heights = defaultdict(lambda: 0)
    for tup in intervals:
        heights[tup[0]] += 1
        heights[tup[1]] -= 1

    retval = []
    height = 0
    left = min(heights.keys())
    for end in sorted(heights.keys()):
        height += heights[end]
        if left is None:
            left = end
        if height == 0:
            retval.append((left, end))
            left = None
    return retval


def save_mmaps(mmaps, buf, total_size):
    nmmaps = unpack_from(buf, 0)
    if nmmaps == 0:
        return
    tuples = ibuf_to_list(buf[sizeof(c_uint64):], 2 * nmmaps, total_size)
    mmaps.extend(tuples)
    pack_to(buf, 0, 0)


def save_preads(preads, buf, idx, total_size):
    npreads = unpack_from(buf, 0)
    if npreads == 0:
        return
    tuples = ibuf_to_list(buf[sizeof(c_uint64):], 2 * npreads, total_size)
    preads[idx] = tuples
    pack_to(buf, 0, 0)


def preprocess_sra(cloud_object: CloudObject):
    logger.info('Preprocessing sra started')

    mapping = get_ncbi_vdb_mapping()
    mapping.dp_sra_size = cloud_object.size
    import tqdm

    mmaps, preads = [], {}
    with (NcbiVdbSharedMemory(mapping.shm_buf, cloud_object.size) as nvsm,
          NcbiVdbSharedMemory(mapping.mmap_buf, MMAP_BUF_SIZE) as mmapsm,
          NcbiVdbSharedMemory(mapping.pread_buf, PREAD_BUF_SIZE) as prsm):
        download(cloud_object, nvsm.buf)
        with temporary_sra() as fpath:
            qc = QuadrupleCursor.from_filepath(fpath)
            total_lines = len(qc)
            save_mmaps(mmaps, mmapsm.buf, cloud_object.size)
            save_preads(preads, prsm.buf, -1, cloud_object.size)
            for i, row in tqdm.tqdm(enumerate(QCRange(qc, total_lines)), total=total_lines):
                save_mmaps(mmaps, mmapsm.buf, cloud_object.size)
                save_preads(preads, prsm.buf, i, cloud_object.size)

    mmaps = sum_intervals(mmaps)
    preads[-1] = sum_intervals(preads[-1])
    return PreprocessingMetadata(
        metadata=io.BytesIO(b'nempty'),
        attributes={
            'total_lines': total_lines,
            'mmaps': mmaps,
            'preads': preads,
        }
    )


@CloudDataFormat(preprocessing_function=preprocess_sra)
class SRA:
    total_lines: int
    index_key: str


class SRASlice(CloudObjectSlice):
    def __init__(self, start, end, mmaps, preads, *args, **kwargs):
        self.start = start
        self.end = end
        self.mmaps = mmaps
        self.preads = preads
        super().__init__(*args, **kwargs)

    def partial_download(self):

        def prefetch(nvsm, mmapsm, prsm):
            nvms_buf_idx = 0
            mmap_buf_idx = 0
            dsize = sizeof(c_uint64)
            pack_to(mmapsm.buf, 0, len(self.mmaps))
            mmap_buf_idx += dsize
            for left, right in self.mmaps:
                length = right - left
                nvsm.buf[nvms_buf_idx:nvms_buf_idx + length] = self.download_range(left, right - 1)
                pack_to(mmapsm.buf, mmap_buf_idx, nvms_buf_idx)
                mmap_buf_idx += dsize
                pack_to(mmapsm.buf, mmap_buf_idx, left)
                mmap_buf_idx += dsize
                pack_to(mmapsm.buf, mmap_buf_idx, right)
                mmap_buf_idx += dsize
                nvms_buf_idx += length

            pread_buf_idx = 0
            pack_to(prsm.buf, 0, len(self.preads))
            pread_buf_idx += dsize
            for left, right in self.preads:
                length = right - left
                nvsm.buf[nvms_buf_idx:nvms_buf_idx + length] = self.download_range(left, right - 1)
                pack_to(prsm.buf, pread_buf_idx, nvms_buf_idx)
                pread_buf_idx += dsize
                pack_to(prsm.buf, pread_buf_idx, left)
                pread_buf_idx += dsize
                pack_to(prsm.buf, pread_buf_idx, right)
                pread_buf_idx += dsize
                nvms_buf_idx += length

        return self.decompress(
            dp_mode=1,
            total_length=self.sum_lengths(self.mmaps) + self.sum_lengths(self.preads),
            prefetch=prefetch
        )


    def full_download(self):

        def prefetch(nvsm, mmapsm, prsm):
            download(self.cloud_object, nvsm.buf)

        return self.decompress(
            dp_mode=0,
            total_length=self.cloud_object.size,
            prefetch=prefetch
        )

    def decompress(self, dp_mode, total_length, prefetch):
        mapping = get_ncbi_vdb_mapping()
        mapping.dp_sra_size = self.cloud_object.size
        mapping.dp_mode = dp_mode
        with (NcbiVdbSharedMemory(mapping.shm_buf, total_length) as nvsm,
              NcbiVdbSharedMemory(mapping.mmap_buf, MMAP_BUF_SIZE) as mmapsm,
              NcbiVdbSharedMemory(mapping.pread_buf, PREAD_BUF_SIZE) as prsm):

            prefetch(nvsm, mmapsm, prsm)

            lines = []
            with temporary_sra() as fpath:
                qc = QuadrupleCursor.from_filepath(fpath)
                for row in QCRange(qc, self.start, self.end):
                    lines.extend(row.lines)
            return lines

    def get(self) -> list[str]:
        return self.partial_download()


    @staticmethod
    def sum_lengths(iterable):
        length = 0
        for left, right in iterable:
            length += right - left

        return length

    def download_range(self, left, right):
        res = self.cloud_object.storage.get_object(
            Bucket=self.cloud_object.path.bucket,
            Key=self.cloud_object.path.key,
            Range=f"bytes={left}-{right}",
        )
        return res['Body'].read()


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


def generate_slices(mmaps, preads, ranges):
    global_preads = preads[-1] + preads[0]
    del preads[-1]
    idx = 0
    slices = []
    local_preads = global_preads.copy()
    prev = None
    for line in sorted(preads.keys()):
        if line >= ranges[idx][1]:
            slices.append(SRASlice(*ranges[idx], mmaps, sum_intervals(local_preads)))
            local_preads = global_preads.copy() + preads[prev]
            idx += 1
        local_preads += preads[line]
        prev = line
    slices.append(SRASlice(*ranges[idx], mmaps, sum_intervals(local_preads)))
    return slices


@PartitioningStrategy(dataformat=SRA)
def partition_chunks_strategy(
    cloud_object: CloudObject,
    num_chunks: int
) -> List[SRASlice]:
    logger.info('SRA partitioning started')
    total_lines = int(cloud_object.get_attribute("total_lines"))
    mmaps = cloud_object.get_attribute("mmaps")
    preads = cloud_object.get_attribute("preads")

    ranges = partition_into_ranges(total_lines, num_chunks)

    return generate_slices(mmaps, preads, ranges)
