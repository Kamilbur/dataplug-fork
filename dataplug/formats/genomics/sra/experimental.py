from __future__ import annotations

import ctypes as C
import io
import logging
import sysconfig
import tempfile
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

from cdlml import get_var

from dataplug.entities import CloudDataFormat, CloudObjectSlice, PartitioningStrategy
from dataplug.preprocessing.metadata import PreprocessingMetadata

from .internals.download import download_range, stream_to_bytes
from .internals.interval import Interval, merge_intervals, partition_dry

if TYPE_CHECKING:
    from dataplug.cloudobject import CloudObject

    from .internals.sra import SRALines

logger = logging.getLogger(__name__)

_types_map = {1: C.c_uint8, 2: C.c_uint16, 4: C.c_uint32, 8: C.c_uint64}
c_off_t = _types_map[sysconfig.get_config_var('SIZEOF_OFF_T')]

MAX_RANGES = 1024
_RangesArray = C.c_uint64 * (1 + MAX_RANGES * 3)

_SYSCALLS = ["open", "close", "fstat", "read", "pread", "mmap", "munmap"]


def _vdb():
    from .internals.vdb import vdb
    return vdb


class FileInfo(C.Structure):
    _fields_ = (
        ('accession', C.c_char_p),
        ('data', C.c_char_p),
        ('size', C.c_size_t),
        ('offset', c_off_t),
    )


class ShimsMapping:
    def __init__(self):
        self._info = get_var(_vdb()['shims'], FileInfo, 'info')

    @property
    def info(self):
        return self._info

    @info.setter
    def info(self, value):
        dst = C.addressof(self._info)
        src = C.addressof(value)
        C.memmove(dst, src, C.sizeof(FileInfo))
        self._keepalive = value


class EnabledMask:
    def __init__(self):
        self._vars = [get_var(_vdb()['shims'], C.c_int, f'enable_{s}') for s in _SYSCALLS]

    @contextmanager
    def enabled_all(self):
        saved = [v.value for v in self._vars]
        for v in self._vars:
            v.value = 1
        try:
            yield
        finally:
            for v, s in zip(self._vars, saved, strict=False):
                v.value = s


_shims_vars = None


def _sv():
    global _shims_vars
    if _shims_vars is None:
        lib = _vdb()['shims']
        _shims_vars = {
            'dp_mode':      get_var(lib, C.c_int, 'dp_mode'),
            'dp_sra_size':  get_var(lib, C.c_size_t, 'dp_sra_size'),
            'pread_ranges': get_var(lib, _RangesArray, 'pread_ranges'),
            'mmap_ranges':  get_var(lib, _RangesArray, 'mmap_ranges'),
        }
    return _shims_vars


def _read_ranges(arr, total_size):
    count = arr[0]
    if count == 0:
        return []
    result = []
    for i in range(count):
        pos = arr[2 * i + 1]
        size = arr[2 * i + 2]
        if pos + size + 4 < total_size:
            size += 4
        result.append((pos, pos + size))
    arr[0] = 0
    return result


def _save_mmaps(mmaps, total_size):
    mmaps.extend(_read_ranges(_sv()['mmap_ranges'], total_size))


def _save_preads(preads, idx, total_size):
    intervals = _read_ranges(_sv()['pread_ranges'], total_size)
    if intervals:
        preads[idx] = intervals


def _fill_ranges_mode1(arr, entries):
    arr[0] = len(entries)
    for i, (buf_off, start, end) in enumerate(entries):
        arr[3 * i + 1] = buf_off
        arr[3 * i + 2] = start
        arr[3 * i + 3] = end


@contextmanager
def temporary_sra(accession):
    with tempfile.TemporaryDirectory() as dpath:
        fpath = Path(dpath) / accession
        fpath.touch()
        yield str(fpath)


def accession(cloud_object: CloudObject) -> str:
    return cloud_object.path.key.split('/')[-1].split('.')[0]


def _download_to_bytes(cloud_object: CloudObject) -> bytes:
    resp = cloud_object.storage.get_object(
        Bucket=cloud_object.path.bucket,
        Key=cloud_object.path.key,
    )
    return stream_to_bytes(resp['Body'])


def _build_range_buffer(cloud_object, mmaps, preads):
    chunks = []
    buf_offset = 0
    mmap_entries = []
    pread_entries = []
    mmaps = merge_intervals(mmaps) if mmaps else []
    preads = merge_intervals(preads) if preads else []
    for left, right in mmaps:
        data = download_range(cloud_object, left, right - 1)
        mmap_entries.append((buf_offset, left, right))
        buf_offset += len(data)
        chunks.append(data)
    for left, right in preads:
        data = download_range(cloud_object, left, right - 1)
        pread_entries.append((buf_offset, left, right))
        buf_offset += len(data)
        chunks.append(data)
    return b''.join(chunks), mmap_entries, pread_entries


def preprocess_sra(cloud_object: CloudObject, step=250):
    from .internals.sra import VColumns
    from .internals.vdb_types import to_char_p

    logger.info('Preprocessing sra started')

    sv = _sv()
    sv['dp_mode'].value = 0
    sv['dp_sra_size'].value = cloud_object.size
    sv['pread_ranges'][0] = 0
    sv['mmap_ranges'][0] = 0

    acc = accession(cloud_object)
    raw = _download_to_bytes(cloud_object)

    shims_mapping = ShimsMapping()
    mask = EnabledMask()
    mmaps, preads = [], {}

    try:
        with mask.enabled_all(), temporary_sra(acc) as fpath:
            shims_mapping.info = FileInfo(
                accession=to_char_p(acc),
                data=to_char_p(raw),
                size=len(raw),
                offset=0,
            )
            with VColumns.from_filepath(fpath) as vcols:
                total_lines = len(vcols)
                _save_mmaps(mmaps, cloud_object.size)
                _save_preads(preads, -1, cloud_object.size)
                for i, row_idx in enumerate(range(1, total_lines + 1, step)):
                    [col.read(row_idx) for col in vcols.columns]
                    _save_mmaps(mmaps, cloud_object.size)
                    _save_preads(preads, i * step, cloud_object.size)
    finally:
        sv['dp_mode'].value = 0

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


@dataclass
class Mapping:
    mmaps: list[Interval] = field(default_factory=list)
    preads: dict[int, list[Interval]] = field(default_factory=dict)

    @classmethod
    def from_cloud_object(cls, cloud_object: CloudObject) -> Mapping:
        mmaps = cloud_object.get_attribute("mmaps")
        preads = cloud_object.get_attribute("preads")
        return cls(mmaps=mmaps, preads=preads)


class SRASlice(CloudObjectSlice):
    def __init__(self, start, end, mmaps, preads, *args, **kwargs):
        self.start = start
        self.end = end
        self.mmaps = mmaps
        self.preads = preads
        super().__init__(*args, **kwargs)

    def partial_download(self, split=False) -> SRALines:
        from .internals.sra import VColumns, _format_spot
        from .internals.vdb_types import to_char_p

        acc = accession(self.cloud_object)
        raw, mmap_entries, pread_entries = _build_range_buffer(
            self.cloud_object, self.mmaps, self.preads
        )
        sv = _sv()
        sv['dp_mode'].value = 1
        sv['dp_sra_size'].value = self.cloud_object.size
        _fill_ranges_mode1(sv['mmap_ranges'], mmap_entries)
        _fill_ranges_mode1(sv['pread_ranges'], pread_entries)

        shims_mapping = ShimsMapping()
        mask = EnabledMask()

        try:
            with mask.enabled_all(), temporary_sra(acc) as fpath:
                shims_mapping.info = FileInfo(
                    accession=to_char_p(acc),
                    data=to_char_p(raw),
                    size=self.cloud_object.size,
                    offset=0,
                )
                with VColumns.from_filepath(fpath) as vcols:
                    for row_idx in range(self.start + 1, self.end + 1):
                        row = [col.read(row_idx) for col in vcols.columns]
                        yield _format_spot(acc, row_idx, *row, split=split)
        finally:
            sv['dp_mode'].value = 0

    def get(self, split=False) -> SRALines:
        return list(self.partial_download(split=split))

    def lazy_get(self, split=False):
        yield from self.partial_download(split=split)

    @staticmethod
    def sum_lengths(iterable):
        return sum(right - left for left, right in iterable)


def generate_slices(mapping: Mapping, ranges: list[Interval]) -> list[SRASlice]:
    mmaps = mapping.mmaps
    preads = {int(key): list(value) for key, value in mapping.preads.items()}
    global_preads = preads.pop(-1, [])
    slices = []
    accumulator = global_preads.copy()

    group_ends = iter(sorted(preads.keys()))
    group_end = next(group_ends, None)

    for range_start, range_end in ranges:
        while group_end is not None and group_end < range_end:
            accumulator.extend(preads[group_end])
            group_end = next(group_ends, None)

        if group_end is not None:
            accumulator.extend(preads[group_end])
            group_end = next(group_ends, None)

        slices.append(SRASlice(
            range_start,
            range_end,
            mmaps,
            merge_intervals(accumulator) if accumulator else [],
        ))

    assert len(slices) == len(ranges)
    return slices


@PartitioningStrategy(dataformat=SRA)
def partition_chunks_strategy(
    cloud_object: CloudObject,
    num_chunks: int,
) -> list[SRASlice]:
    logger.info('SRA partitioning started')
    total_lines = int(cloud_object.get_attribute("total_lines"))
    mapping = Mapping.from_cloud_object(cloud_object)
    ranges = partition_dry(total_lines, num_chunks)
    return generate_slices(mapping, ranges)
