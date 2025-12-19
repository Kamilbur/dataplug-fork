from __future__ import annotations

from contextlib import contextmanager
from ctypes import (
    c_uint64,
    sizeof,
)
from dataclasses import dataclass, field
import io
import logging
import os
from pathlib import Path
import tempfile
from typing import TYPE_CHECKING

from dataplug.entities import CloudDataFormat, CloudObjectSlice, PartitioningStrategy
from dataplug.preprocessing.metadata import PreprocessingMetadata

from .internals._loader import get_table
from .internals.download import download, download_range
from .internals.experimental.cmaps import (
    NcbiVdbSharedMemory,
    get_ncbi_vdb_mapping,
    pack_to,
    save_mmaps,
    save_preads,
)
from .internals.experimental.interval import Interval, merge_intervals, partition_dry
from .internals.sra import SRAFile
from .output import FastqLineQuadruple

logger = logging.getLogger(__name__)

PREAD_BUF_SIZE = 100 * sizeof(c_uint64)
MMAP_BUF_SIZE = 100 * sizeof(c_uint64)

libvdb_python_path = 'LIBVDB_PYTHON_DIR'
libvdb_path = os.environ.get('NCBI_VDB_SO_PATH', '')

if not libvdb_path:
    logger.warning(
        'For complete functionality of SRA backend provide path'
        'to ncbi-vdb shared object in NCBI_VDB_SO_PATH environment variable.'
    )


if TYPE_CHECKING:

    from dataplug.cloudobject import CloudObject

    from .output import SRALines


# For some reason, shared buffer does not work with default implementation,
# so we create simplified one here
@dataclass
class QuadrupleCursor:
    columns: dict

    _column_names = (
        "READ",
        "(INSDC:quality:text:phred_33)QUALITY",
        "NAME",
        "READ_START",
        "READ_LEN",
    )

    @classmethod
    def from_filepath(cls, filepath):
        table = get_table(filepath)

        columns = table.CreateCursor().OpenColumns(list(cls._column_names))

        return cls(columns)

    def row(self, idx):
        return [col.Read(idx) for col in self.columns.values()]

    def __len__(self):
        return self.columns[self._column_names[0]].row_range()[1]


def qc_range(cursor, start, stop, step):
    for i in range(start + 1, stop + 1, step):
        yield cursor.row(i)


@contextmanager
def temporary_sra(accession):
    with tempfile.TemporaryDirectory() as dpath:
        fpath = Path(dpath) / accession
        with Path(fpath).open('w+', encoding='utf-8') as fp:
            fp.write('s3re')
            fp.flush()
            yield str(fpath)


def accession(cloud_object: CloudObject) -> str:
    return cloud_object.path.key.split('/')[-1].split('.')[0]


@dataclass
class Mapping:
    mmaps: list[Interval] = field(default_factory=list)
    preads: dict[int, list[Interval]] = field(default_factory=dict)

    @classmethod
    def from_cloud_object(cls, cloud_object: CloudObject) -> Mapping:
        mmaps = cloud_object.get_attribute("mmaps")
        preads = cloud_object.get_attribute("preads")
        return cls(mmaps=mmaps, preads=preads)

    def save(self, mmaps, preads, idx, total_size):
        save_mmaps(mmaps, mmaps, total_size)
        save_preads(preads, preads, idx, total_size)


def preprocess_sra(cloud_object: CloudObject, step=250):
    logger.info('Preprocessing sra started')

    ncbi_mapping = get_ncbi_vdb_mapping(libvdb_path)
    ncbi_mapping.dp_sra_size = cloud_object.size
    ncbi_mapping.dp_mode = 0

    mapping = Mapping()

    acc = accession(cloud_object)
    mmaps, preads = [], {}
    with (NcbiVdbSharedMemory(ncbi_mapping.shm_buf, cloud_object.size) as nvsm,
          NcbiVdbSharedMemory(ncbi_mapping.mmap_buf, MMAP_BUF_SIZE) as mmapsm,
          NcbiVdbSharedMemory(ncbi_mapping.pread_buf, PREAD_BUF_SIZE) as prsm,
                temporary_sra(acc) as fpath):
        download(cloud_object, nvsm.buf)
        qc = QuadrupleCursor.from_filepath(fpath)
        total_lines = len(qc)
        # mapping.save(mmaps, preads, -1, cloud_object.size)
        save_mmaps(mmaps, mmapsm.buf, cloud_object.size)
        save_preads(preads, prsm.buf, -1, cloud_object.size)
        for i, _ in enumerate(qc_range(qc, 0, total_lines, step=step)):
            # mapping.save(mmaps, preads, i * step, cloud_object.size)
            save_mmaps(mmaps, mmapsm.buf, cloud_object.size)
            save_preads(preads, prsm.buf, i * step, cloud_object.size)

    # mmaps = merge_intervals(mmaps)

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


def prefetch(mmaps, preads, cloud_object, nvsm, mmapsm, prsm):
    nvms_buf_idx = 0
    mmap_buf_idx = 0
    dsize = sizeof(c_uint64)
    pack_to(mmapsm.buf, 0, len(mmaps))
    mmap_buf_idx += dsize
    for left, right in mmaps:
        length = right - left
        nvsm.buf[nvms_buf_idx:nvms_buf_idx + length] = download_range(cloud_object, left, right - 1)
        pack_to(mmapsm.buf, mmap_buf_idx, nvms_buf_idx)
        mmap_buf_idx += dsize
        pack_to(mmapsm.buf, mmap_buf_idx, left)
        mmap_buf_idx += dsize
        pack_to(mmapsm.buf, mmap_buf_idx, right)
        mmap_buf_idx += dsize
        nvms_buf_idx += length

    pread_buf_idx = 0
    pack_to(prsm.buf, 0, len(preads))
    pread_buf_idx += dsize
    for left, right in preads:
        length = right - left
        nvsm.buf[nvms_buf_idx:nvms_buf_idx + length] = download_range(cloud_object, left, right - 1)
        pack_to(prsm.buf, pread_buf_idx, nvms_buf_idx)
        pread_buf_idx += dsize
        pack_to(prsm.buf, pread_buf_idx, left)
        pread_buf_idx += dsize
        pack_to(prsm.buf, pread_buf_idx, right)
        pread_buf_idx += dsize
        nvms_buf_idx += length


class SRASlice(CloudObjectSlice):
    def __init__(self, start, end, mmaps, preads, *args, **kwargs):
        self.start = start
        self.end = end
        self.mmaps = mmaps
        self.preads = preads
        super().__init__(*args, **kwargs)

    def partial_download(self):
        yield from self.decompress(
            dp_mode=1,
            total_length=self.sum_lengths(self.mmaps) + self.sum_lengths(self.preads),
            prefetch=prefetch
        )

    def decompress(self, dp_mode, total_length, prefetch):
        mapping = get_ncbi_vdb_mapping()
        mapping.dp_sra_size = self.cloud_object.size
        mapping.dp_mode = dp_mode
        with (NcbiVdbSharedMemory(mapping.shm_buf, total_length) as nvsm,
              NcbiVdbSharedMemory(mapping.mmap_buf, MMAP_BUF_SIZE) as mmapsm,
              NcbiVdbSharedMemory(mapping.pread_buf, PREAD_BUF_SIZE) as prsm,
                temporary_sra(accession(self.cloud_object)) as fpath):
            prefetch(self.mmaps, self.preads, self.cloud_object, nvsm, mmapsm, prsm)
            sra_file = SRAFile(fpath, paired=False)
            for spot in sra_file.range(self.start, self.end):
                reads_repr = tuple(str(FastqLineQuadruple(read)) for read in spot.reads())
                write_unpaired = False
                if not write_unpaired:
                    yield reads_repr
                else:
                    yield reads_repr

    def get(self) -> SRALines:
        return list(self.partial_download())

    def lazy_get(self):
        yield from self.partial_download()

    @staticmethod
    def sum_lengths(iterable):
        length = 0
        for left, right in iterable:
            length += right - left

        return length


# def generate_slices(mapping: Mapping, ranges: list[Interval]) -> list[SRASlice]:
#     preads, mmaps = mapping.preads, mapping.mmaps
#     global_preads = preads[-1] + preads[0]
#     del preads[-1]
#     idx = 0
#     slices = []
#     local_preads = global_preads.copy()
#     prev = -1

#     for line in sorted(preads.keys()):
#         if line >= ranges[idx][1]:
#             local_preads += preads[line]
#             while line >= ranges[idx][1]:
#                 slices.append(SRASlice(*ranges[idx], mmaps, merge_intervals(local_preads)))
#                 idx += 1
#             local_preads = global_preads.copy() + preads[prev]

#         local_preads += preads[line]
#         prev = line

#     while len(slices) < len(ranges):
#         slices.append(SRASlice(*ranges[idx], mmaps, merge_intervals(local_preads)))
#         idx += 1

#     assert len(slices) == len(ranges)
#     return slices


def generate_slices(mapping: Mapping, ranges: list[Interval]) -> list[SRASlice]:
    preads, mmaps = mapping.preads, mapping.mmaps
    global_preads = preads[-1]# + preads[0]
    del preads[-1]
    i = 0
    j = 0
    slices = []
    accumulator = global_preads.copy()

    group_ends = sorted(preads.keys())
    while i < len(group_ends) and j < len(ranges):
        g = group_ends[i]
        p = ranges[j][1]
        accumulator += preads[g]

        if g < p:
            i += 1
        else:
            slices.append(SRASlice(*ranges[j], mmaps, merge_intervals(accumulator)))
            j += 1
            if j < len(ranges):
                accumulator = global_preads.copy() #+ preads[G[i - 1]]

    while j < len(ranges):
        slices.append(SRASlice(*ranges[j], mmaps, merge_intervals(accumulator)))
        j += 1
    assert len(slices) == len(ranges)
    return slices


@PartitioningStrategy(dataformat=SRA)
def partition_chunks_strategy(
    cloud_object: CloudObject,
    num_chunks: int
) -> list[SRASlice]:
    logger.info('SRA partitioning started')

    total_lines = int(cloud_object.get_attribute("total_lines"))
    mapping = Mapping.from_cloud_object(cloud_object)

    ranges = partition_dry(total_lines, num_chunks)
    return generate_slices(mapping, ranges)
