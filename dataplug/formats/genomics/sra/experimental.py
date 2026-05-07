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

from .internals.download import download_range, parallel_download_to_buffer, stream_to_bytes
from .internals.interval import Interval, merge_intervals, partition_dry

if TYPE_CHECKING:
    from dataplug.cloudobject import CloudObject

    from .internals.sra import SRALines

logger = logging.getLogger(__name__)

_types_map = {1: C.c_uint8, 2: C.c_uint16, 4: C.c_uint32, 8: C.c_uint64}
c_off_t = _types_map[sysconfig.get_config_var("SIZEOF_OFF_T")]

MAX_RANGES = 1024
_RangesArray = C.c_uint64 * (1 + MAX_RANGES * 3)

_SYSCALLS = ["open", "close", "fstat", "read", "pread", "mmap", "munmap", "lseek"]


def _vdb():
    from .internals.vdb import vdb

    return vdb


def _reset_vdb():
    global _shims_vars
    from .internals.vdb import reset_vdb

    _shims_vars = None
    return reset_vdb()


class FileInfo(C.Structure):
    _fields_ = (
        ("accession", C.c_char_p),
        ("data", C.c_char_p),
        ("size", C.c_size_t),
        ("offset", c_off_t),
    )


def to_bytes(value, size: int | None = None) -> bytes:
    if isinstance(value, bytes):
        return value if size is None else value[:size]
    if isinstance(value, str):
        data = value.encode()
        return data if size is None else data[:size]
    if isinstance(value, C.c_char_p):
        if size is None:
            return value.value or b""
        return C.string_at(value, size)
    if size is None:
        return bytes(value)
    return C.string_at(value, size)


class ShimsMapping:
    def __init__(self):
        self._lib = _vdb()["shims"]
        self._set_info = self._lib.dp_set_info_buffer
        self._set_info.argtypes = [C.c_char_p, C.c_char_p, C.c_size_t, C.c_size_t, c_off_t]
        self._set_info.restype = None
        self._set_info_zc = self._lib.dp_set_info_zerocopy
        self._set_info_zc.argtypes = [C.c_char_p, C.c_void_p, C.c_size_t, C.c_size_t, c_off_t]
        self._set_info_zc.restype = None
        self._clear_info = self._lib.dp_clear_info
        self._clear_info.argtypes = []
        self._clear_info.restype = None

    @property
    def info(self):
        return get_var(self._lib, FileInfo, "info")

    @info.setter
    def info(self, value):
        self._keepalive = value
        self.set_info(value.accession, value.data, value.size, value.offset)

    def set_info(self, accession, data, size, offset=0):
        self._keepalive = (accession, data)
        data_bytes = to_bytes(data)
        self._set_info(to_bytes(accession), data_bytes, len(data_bytes), size, offset)

    def set_info_zerocopy(self, accession, data, size, offset=0):
        """Pass `data` to the C shim without copying. `data` must outlive the
        next clear() call; we keep a strong reference here.
        """
        if hasattr(self._lib, "_server"):
            self.set_info(accession, data, size, offset)
            return

        self._keepalive = (accession, data)
        if isinstance(data, (bytes, bytearray, memoryview)):
            buf = data if isinstance(data, bytearray) else bytearray(data)
            ptr = (C.c_char * len(buf)).from_buffer(buf)
            self._keepalive = (accession, buf, ptr)
            addr = C.addressof(ptr)
        else:
            addr = data
            buf = data
        self._set_info_zc(to_bytes(accession), addr, len(buf), size, offset)


    def clear(self):
        self._clear_info()
        self._keepalive = None


class EnabledMask:
    def __init__(self):
        self._vars = [get_var(_vdb()["shims"], C.c_int, f"enable_{s}") for s in _SYSCALLS]

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
        lib = _vdb()["shims"]
        _shims_vars = {
            "dp_mode": get_var(lib, C.c_int, "dp_mode"),
            "dp_sra_size": get_var(lib, C.c_size_t, "dp_sra_size"),
            "pread_ranges": get_var(lib, _RangesArray, "pread_ranges"),
            "mmap_ranges": get_var(lib, _RangesArray, "mmap_ranges"),
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
    mmaps.extend(_read_ranges(_sv()["mmap_ranges"], total_size))


def _save_preads(preads, idx, total_size):
    intervals = _read_ranges(_sv()["pread_ranges"], total_size)
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
    return cloud_object.path.key.split("/")[-1].split(".")[0]


def _download_to_bytes(cloud_object: CloudObject) -> bytes:
    resp = cloud_object.storage.get_object(
        Bucket=cloud_object.path.bucket,
        Key=cloud_object.path.key,
    )
    return stream_to_bytes(resp["Body"])


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
    return b"".join(chunks), mmap_entries, pread_entries


def _walk_record_ranges_python(vcols, total_lines: int, step: int, total_size: int,
                               mmaps: list, preads: dict) -> None:
    """Walk every `step` rows, triggering VDB column reads (without
    materializing data) and capturing the mmap/pread ranges that the shim
    records for each row group.
    """
    import numpy as np

    cur = vcols.columns[0].cur
    col_idxs = tuple(C.c_int(col.idx) for col in vcols.columns)
    elem_bits = C.c_int()
    data = C.c_void_p()
    row_len = C.c_int()
    eb_ref = C.byref(elem_bits)
    d_ref = C.byref(data)
    rl_ref = C.byref(row_len)
    rid = C.c_longlong(0)
    cell = _vdb().VCursorCellDataDirect

    sv = _sv()
    mmap_arr = sv["mmap_ranges"]
    pread_arr = sv["pread_ranges"]
    mmap_view = np.frombuffer(mmap_arr, dtype=np.uint64)
    pread_view = np.frombuffer(pread_arr, dtype=np.uint64)

    def drain(arr_view) -> list:
        count = int(arr_view[0])
        if count == 0:
            return []
        body = arr_view[1 : 1 + 2 * count].reshape(-1, 2)
        positions = body[:, 0]
        sizes = body[:, 1].copy()
        room = positions + sizes + 4 < total_size
        sizes[room] += 4
        ends = positions + sizes
        result = list(zip(positions.tolist(), ends.tolist()))
        arr_view[0] = 0
        return result

    mmaps.extend(drain(mmap_view))
    initial_pread = drain(pread_view)
    if initial_pread:
        preads[-1] = initial_pread

    i = 0
    row_idx = 1
    while row_idx <= total_lines:
        rid.value = row_idx
        for ci in col_idxs:
            cell(cur, rid, ci, eb_ref, d_ref, None, rl_ref)
        m = drain(mmap_view)
        if m:
            mmaps.extend(m)
        p = drain(pread_view)
        if p:
            preads[i * step] = p
        i += 1
        row_idx += step


def _walk_record_ranges_c(vcols, total_lines: int, step: int, total_size: int,
                          mmaps: list, preads: dict) -> bool:
    import numpy as np
    from .internals.cdlml_compat import cast, pointer_for, string_at

    lib = _vdb()["shims"]
    try:
        set_cell = lib.dp_set_vcursor_cell_func
        walk = lib.dp_walk_columns
        clear = lib.dp_clear_walk_ranges
        pread_count = lib.dp_walk_pread_count
        pread_data = lib.dp_walk_pread_data
        mmap_count = lib.dp_walk_mmap_count
        mmap_data = lib.dp_walk_mmap_data
    except AttributeError:
        return False

    set_cell.argtypes = [C.c_void_p]
    set_cell.restype = None
    set_cell(C.cast(_vdb().VCursorCellDataDirect, C.c_void_p))
    walk.argtypes = [
        C.c_void_p,
        C.POINTER(C.c_int),
        C.c_size_t,
        C.c_longlong,
        C.c_longlong,
        C.c_uint64,
    ]
    walk.restype = C.c_int
    clear.argtypes = []
    clear.restype = None
    pread_count.argtypes = []
    pread_count.restype = C.c_size_t
    pread_data.argtypes = []
    pread_data.restype = C.POINTER(C.c_uint64)
    mmap_count.argtypes = []
    mmap_count.restype = C.c_size_t
    mmap_data.argtypes = []
    mmap_data.restype = C.POINTER(C.c_uint64)

    col_idxs = (C.c_int * len(vcols.columns))(*(col.idx for col in vcols.columns))
    col_idxs_arg = pointer_for(col_idxs, cdll=lib) if hasattr(lib, "_server") else col_idxs
    cur_arg = vcols.cur.value if isinstance(vcols.cur, C.c_void_p) else vcols.cur
    rc = walk(cur_arg, col_idxs_arg, len(col_idxs), total_lines, step, total_size)
    if rc == -10:
        clear()
        return False
    if rc != 0:
        clear()
        raise RuntimeError(f"dp_walk_columns failed with rc={rc}")

    uint64_max = (1 << 64) - 1
    try:
        mc = mmap_count()
        if mc:
            data = mmap_data()
            if hasattr(lib, "_server"):
                data = string_at(cast(data, C.POINTER(C.c_uint64), cdll=lib), mc * 2 * C.sizeof(C.c_uint64))
                arr = np.frombuffer(data, dtype=np.uint64).reshape(-1, 2)
            else:
                arr = np.ctypeslib.as_array(data, shape=(mc * 2,)).reshape(-1, 2)
            mmaps.extend((int(row[0]), int(row[1])) for row in arr)

        pc = pread_count()
        if pc:
            data = pread_data()
            if hasattr(lib, "_server"):
                data = string_at(cast(data, C.POINTER(C.c_uint64), cdll=lib), pc * 3 * C.sizeof(C.c_uint64))
                arr = np.frombuffer(data, dtype=np.uint64).reshape(-1, 3)
            else:
                arr = np.ctypeslib.as_array(data, shape=(pc * 3,)).reshape(-1, 3)
            for row in arr:
                key = int(row[0])
                idx = -1 if key == uint64_max else int(key)
                preads.setdefault(idx, []).append((int(row[1]), int(row[2])))
    finally:
        clear()
    return True


def _walk_record_ranges(vcols, total_lines: int, step: int, total_size: int,
                         mmaps: list, preads: dict) -> None:
    if _walk_record_ranges_c(vcols, total_lines, step, total_size, mmaps, preads):
        return
    _walk_record_ranges_python(vcols, total_lines, step, total_size, mmaps, preads)


def preprocess_sra(cloud_object: CloudObject, step=250):
    from .internals.sra import VColumns

    logger.info("Preprocessing sra started")

    _reset_vdb()
    sv = _sv()
    sv["dp_mode"].value = 0
    sv["dp_sra_size"].value = cloud_object.size
    sv["pread_ranges"][0] = 0
    sv["mmap_ranges"][0] = 0

    acc = accession(cloud_object)
    raw = parallel_download_to_buffer(cloud_object)

    shims_mapping = ShimsMapping()
    mask = EnabledMask()
    mmaps, preads = [], {}

    try:
        with mask.enabled_all(), temporary_sra(acc) as fpath:
            shims_mapping.set_info_zerocopy(acc, raw, len(raw), 0)
            with VColumns.from_filepath(fpath) as vcols:
                total_lines = len(vcols)
                _walk_record_ranges(
                    vcols, total_lines, step, cloud_object.size, mmaps, preads
                )
    finally:
        sv["dp_mode"].value = 0
        shims_mapping.clear()
        del raw

    return PreprocessingMetadata(
        metadata=io.BytesIO(b"nempty"),
        attributes={
            "total_lines": total_lines,
            "mmaps": mmaps,
            "preads": preads,
        },
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

        acc = accession(self.cloud_object)
        raw, mmap_entries, pread_entries = _build_range_buffer(
            self.cloud_object, self.mmaps, self.preads
        )
        sv = _sv()
        sv["dp_mode"].value = 1
        sv["dp_sra_size"].value = self.cloud_object.size
        _fill_ranges_mode1(sv["mmap_ranges"], mmap_entries)
        _fill_ranges_mode1(sv["pread_ranges"], pread_entries)

        shims_mapping = ShimsMapping()
        mask = EnabledMask()

        try:
            with mask.enabled_all(), temporary_sra(acc) as fpath:
                shims_mapping.set_info(acc, raw, self.cloud_object.size, 0)
                with VColumns.from_filepath(fpath) as vcols:
                    for row_idx in range(self.start + 1, self.end + 1):
                        row = [col.read(row_idx) for col in vcols.columns]
                        yield _format_spot(acc, row_idx, *row, split=split)
        finally:
            sv["dp_mode"].value = 0
            shims_mapping.clear()

    def get(self, split=False) -> SRALines:
        return list(self.partial_download(split=split))

    def lazy_get(self, split=False):
        yield from self.partial_download(split=split)

    def to_file_obj(self, file_obj, split=False):
        if split or not hasattr(file_obj, "fileno"):
            items = 0
            output_bytes = 0
            for reads in self.partial_download(split=split):
                for record in reads:
                    written = file_obj.write(record.encode("utf-8"))
                    output_bytes += written if isinstance(written, int) else len(record)
                    items += 1
            return items, output_bytes

        from .internals.sra import VColumns

        acc = accession(self.cloud_object)
        raw, mmap_entries, pread_entries = _build_range_buffer(
            self.cloud_object, self.mmaps, self.preads
        )
        sv = _sv()
        sv["dp_mode"].value = 1
        sv["dp_sra_size"].value = self.cloud_object.size
        _fill_ranges_mode1(sv["mmap_ranges"], mmap_entries)
        _fill_ranges_mode1(sv["pread_ranges"], pread_entries)

        shims_mapping = ShimsMapping()
        mask = EnabledMask()
        items = C.c_uint64()
        output_bytes = C.c_uint64()

        try:
            with mask.enabled_all(), temporary_sra(acc) as fpath:
                shims_mapping.set_info(acc, raw, self.cloud_object.size, 0)
                with VColumns.from_filepath(fpath) as vcols:
                    lib = _vdb()["shims"]
                    set_cell = lib.dp_set_vcursor_cell_func
                    set_cell.argtypes = [C.c_void_p]
                    set_cell.restype = None
                    set_cell(C.cast(_vdb().VCursorCellDataDirect, C.c_void_p))
                    write_range = lib.dp_write_fastq_range
                    write_range.argtypes = [
                        C.c_void_p,
                        C.c_int,
                        C.c_int,
                        C.c_int,
                        C.c_longlong,
                        C.c_longlong,
                        C.c_char_p,
                        C.c_int,
                        C.POINTER(C.c_uint64),
                        C.POINTER(C.c_uint64),
                    ]
                    write_range.restype = C.c_int
                    cols = vcols.columns
                    rc = write_range(
                        vcols.cur,
                        cols[0].idx,
                        cols[1].idx,
                        cols[2].idx,
                        self.start + 1,
                        self.end,
                        acc.encode(),
                        file_obj.fileno(),
                        C.byref(items),
                        C.byref(output_bytes),
                    )
                    if rc != 0:
                        raise RuntimeError(f"dp_write_fastq_range failed with rc={rc}")
        finally:
            sv["dp_mode"].value = 0
            shims_mapping.clear()
        return int(items.value), int(output_bytes.value)

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

        slices.append(
            SRASlice(
                range_start,
                range_end,
                mmaps,
                merge_intervals(accumulator) if accumulator else [],
            )
        )

    assert len(slices) == len(ranges)
    return slices


@PartitioningStrategy(dataformat=SRA)
def partition_chunks_strategy(
    cloud_object: CloudObject,
    num_chunks: int,
) -> list[SRASlice]:
    logger.info("SRA partitioning started")
    total_lines = int(cloud_object.get_attribute("total_lines"))
    mapping = Mapping.from_cloud_object(cloud_object)
    ranges = partition_dry(total_lines, num_chunks)
    return generate_slices(mapping, ranges)
