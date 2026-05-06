from __future__ import annotations

import ctypes as C
from contextlib import contextmanager
from dataclasses import dataclass, field
import io
import logging
from pathlib import Path
import tempfile
from typing import TYPE_CHECKING

from cdlml import get_var

from dataplug.entities import CloudDataFormat, CloudObjectSlice, PartitioningStrategy
from dataplug.preprocessing.metadata import PreprocessingMetadata

from .internals.download import download_range, stream_to_bytes
from .internals.interval import Interval, merge_intervals, partition_dry

if TYPE_CHECKING:
    from dataplug.cloudobject import CloudObject


logger = logging.getLogger(__name__)

MAX_RANGES = 65536
_RangesArray = C.c_uint64 * (1 + MAX_RANGES * 3)
_SYSCALLS = ["open", "open64", "close", "fstat", "fstat64", "read", "lseek", "lseek64"]


def _shims():
    from .internals.genie import genie
    return genie["shims"]


class FileInfo(C.Structure):
    _fields_ = (
        ("filename", C.c_char_p),
        ("data", C.c_char_p),
        ("size", C.c_size_t),
    )


class ShimsMapping:
    def __init__(self):
        self._info = get_var(_shims(), FileInfo, "info")

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
        self._vars = [get_var(_shims(), C.c_int, f"enable_{syscall}") for syscall in _SYSCALLS]

    @contextmanager
    def enabled_all(self):
        saved = [v.value for v in self._vars]
        for var in self._vars:
            var.value = 1
        try:
            yield
        finally:
            for var, value in zip(self._vars, saved):
                var.value = value


_shims_vars = None


def _sv():
    global _shims_vars
    if _shims_vars is None:
        lib = _shims()
        _shims_vars = {
            "dp_mode": get_var(lib, C.c_int, "dp_mode"),
            "read_ranges": get_var(lib, _RangesArray, "read_ranges"),
        }
    return _shims_vars


def _read_ranges(arr):
    count = arr[0]
    result = []
    for idx in range(count):
        left = arr[3 * idx + 2]
        right = arr[3 * idx + 3]
        if right > left:
            result.append((left, right))
    arr[0] = 0
    return result


def _save_reads(reads: dict[int, list[Interval]], idx: int) -> None:
    intervals = _read_ranges(_sv()["read_ranges"])
    if intervals:
        reads[idx] = merge_intervals(intervals)


def _fill_ranges_mode1(arr, entries) -> None:
    arr[0] = len(entries)
    for idx, (buf_off, start, end) in enumerate(entries):
        arr[3 * idx + 1] = buf_off
        arr[3 * idx + 2] = start
        arr[3 * idx + 3] = end


@contextmanager
def temporary_mgb(filename: str, size: int):
    with tempfile.TemporaryDirectory() as dpath:
        fpath = Path(dpath) / filename
        with fpath.open("wb") as f:
            f.truncate(size)
        yield str(fpath), dpath


@contextmanager
def temporary_mgb_bytes(filename: str, data: bytes):
    with tempfile.TemporaryDirectory() as dpath:
        fpath = Path(dpath) / filename
        fpath.write_bytes(data)
        yield str(fpath), dpath


def _write_ranges(path: str, data: bytes, entries) -> None:
    with open(path, "r+b") as f:
        for buf_off, start, end in entries:
            f.seek(start)
            f.write(data[buf_off:buf_off + (end - start)])


def filename(cloud_object: CloudObject) -> str:
    name = cloud_object.path.key.split("/")[-1]
    return name if name.endswith(".mgb") else f"{name}.mgb"


def _download_to_bytes(cloud_object: CloudObject) -> bytes:
    response = cloud_object.storage.get_object(
        Bucket=cloud_object.path.bucket,
        Key=cloud_object.path.key,
    )
    return stream_to_bytes(response["Body"])


def _build_range_buffer(cloud_object: CloudObject, intervals: list[Interval]):
    chunks = []
    entries = []
    offset = 0
    for left, right in merge_intervals(intervals):
        data = download_range(cloud_object, left, right - 1)
        entries.append((offset, left, right))
        chunks.append(data)
        offset += len(data)
    return b"".join(chunks), entries


def _fastq_records(data: bytes):
    text = data.decode("utf-8")
    lines = text.rstrip("\n").split("\n")
    records = []
    for idx in range(0, len(lines), 4):
        record_lines = lines[idx:idx + 4]
        if len(record_lines) == 4:
            if record_lines[2] == "+":
                record_lines[2] = "+" + record_lines[0][1:]
            records.append(record_lines)

    for record_lines in sorted(records, key=_fastq_record_sort_key):
        yield ("\n".join(record_lines) + "\n",)


def _fastq_record_sort_key(record_lines: list[str]):
    read_name = record_lines[0].split()[0]
    try:
        return (0, int(read_name.rsplit(".", 1)[1]))
    except (IndexError, ValueError):
        return (1, read_name)


def _ensure_read_range(reads: dict[int, list[Interval]], idx: int, size: int) -> None:
    if idx not in reads:
        reads[idx] = [(0, size)]


def preprocess_mgb(cloud_object: CloudObject, threads: int = 1):
    from .internals.genie import access_unit_count, decompress_access_unit

    logger.info("Preprocessing mgb started")
    sv = _sv()
    sv["dp_mode"].value = 0
    sv["read_ranges"][0] = 0

    file_name = filename(cloud_object)
    raw = _download_to_bytes(cloud_object)
    raw_buffer = C.create_string_buffer(raw)
    raw_c = C.cast(raw_buffer, C.c_char_p)
    name_c = C.c_char_p(file_name.encode())

    shims_mapping = ShimsMapping()
    mask = EnabledMask()
    reads: dict[int, list[Interval]] = {}

    try:
        with mask.enabled_all(), temporary_mgb_bytes(file_name, raw) as (fpath, workdir):
            shims_mapping.info = FileInfo(name_c, raw_c, len(raw))
            num_access_units = access_unit_count(fpath)
            _save_reads(reads, -1)
            for access_unit_id in range(num_access_units):
                decompress_access_unit(
                    fpath,
                    access_unit_id,
                    working_dir=workdir,
                    threads=threads,
                )
                _save_reads(reads, access_unit_id)
                _ensure_read_range(reads, access_unit_id, len(raw))
    finally:
        sv["dp_mode"].value = 0

    return PreprocessingMetadata(
        metadata=io.BytesIO(b"nempty"),
        attributes={
            "access_units": num_access_units,
            "reads": reads,
        },
    )


@CloudDataFormat(preprocessing_function=preprocess_mgb)
class MGB:
    access_units: int
    reads: dict[int, list[Interval]]


@dataclass
class Mapping:
    reads: dict[int, list[Interval]] = field(default_factory=dict)

    @classmethod
    def from_cloud_object(cls, cloud_object: CloudObject) -> Mapping:
        return cls(reads=cloud_object.get_attribute("reads"))


class MGBSlice(CloudObjectSlice):
    def __init__(self, start: int, end: int, reads: list[Interval], *args, **kwargs):
        self.start = start
        self.end = end
        self.reads = reads
        super().__init__(*args, **kwargs)

    def partial_download(self, threads: int = 1):
        from .internals.genie import decompress_access_unit

        file_name = filename(self.cloud_object)
        raw, entries = _build_range_buffer(self.cloud_object, self.reads)
        with temporary_mgb(file_name, self.cloud_object.size) as (fpath, workdir):
            _write_ranges(fpath, raw, entries)
            for access_unit_id in range(self.start, self.end):
                data = decompress_access_unit(
                    fpath,
                    access_unit_id,
                    working_dir=workdir,
                    threads=threads,
                )
                yield from _fastq_records(data)

    def get(self, threads: int = 1):
        return list(self.partial_download(threads=threads))

    def lazy_get(self, threads: int = 1):
        yield from self.partial_download(threads=threads)


def generate_slices(mapping: Mapping, ranges: list[Interval]) -> list[MGBSlice]:
    reads = {int(key): list(value) for key, value in mapping.reads.items()}
    global_reads = reads.pop(-1, [])
    slices = []
    for start, end in ranges:
        intervals = list(global_reads)
        for access_unit_id in range(start, end):
            intervals.extend(reads.get(access_unit_id, []))
        slices.append(MGBSlice(start, end, merge_intervals(intervals)))
    return slices


@PartitioningStrategy(dataformat=MGB)
def partition_chunks_strategy(cloud_object: CloudObject, num_chunks: int) -> list[MGBSlice]:
    logger.info("MGB partitioning started")
    total = int(cloud_object.get_attribute("access_units"))
    return generate_slices(Mapping.from_cloud_object(cloud_object), partition_dry(total, num_chunks))
