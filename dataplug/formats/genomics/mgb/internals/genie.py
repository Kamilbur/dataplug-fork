from __future__ import annotations

import ctypes as C
import os
from pathlib import Path

from cdlml import Preload, PreloadedCDLL, byref, cast, string_at

_HERE = Path(__file__).parent
_REPO_ROOT = _HERE.parents[4]

GENIE_SHARED_SUCCESS = 0
GENIE_LOG_DEBUG = 0
GENIE_LOG_INFO = 1
GENIE_LOG_WARNING = 2
GENIE_LOG_ERROR = 3


def _default_genie_path() -> Path:
    env_path = os.environ.get("GENIE_SO_PATH")
    if env_path:
        return Path(env_path)

    local = _HERE / "libgenie.so"
    if local.exists():
        return local

    return local


preloads = [
    Preload(name=b"libmgbshims.so", path=str(_HERE / "libmgbshims.so").encode(), alias="shims"),
]

genie = PreloadedCDLL(str(_default_genie_path()), preloads=preloads)

genie.GenieSharedStrerror.argtypes = [C.c_uint8]
genie.GenieSharedStrerror.restype = C.c_char_p

genie.GenieSetLogSeverity.argtypes = [C.c_uint8]
genie.GenieSetLogSeverity.restype = C.c_uint8

genie.GenieGetAccessUnitCount.argtypes = [C.c_char_p, C.POINTER(C.c_uint64)]
genie.GenieGetAccessUnitCount.restype = C.c_uint8

genie.GenieDecompressAccessUnit.argtypes = [
    C.c_char_p,
    C.c_uint64,
    C.c_char_p,
    C.c_char_p,
    C.c_char_p,
    C.c_uint64,
]
genie.GenieDecompressAccessUnit.restype = C.c_uint8

genie.GenieDecompressAccessUnitToFastq.argtypes = [
    C.c_char_p,
    C.c_uint64,
    C.c_char_p,
    C.c_char_p,
    C.c_uint64,
    C.POINTER(C.c_void_p),
    C.POINTER(C.c_uint64),
]
genie.GenieDecompressAccessUnitToFastq.restype = C.c_uint8

genie.GenieFree.argtypes = [C.c_void_p]
genie.GenieFree.restype = None


class GenieError(RuntimeError):
    def __init__(self, code: int) -> None:
        raw = genie.GenieSharedStrerror(code)
        msg = raw.decode() if raw else "unknown error"
        super().__init__(f"GENIE error {code}: {msg}")
        self.code = code


def _path(path: str | Path | None) -> bytes | None:
    if path is None:
        return None
    return str(path).encode()


def _optional_path(path: str | Path | None) -> bytes:
    if path is None:
        return b""
    return str(path).encode()


def _check(code: int) -> None:
    if code != GENIE_SHARED_SUCCESS:
        raise GenieError(code)


def set_log_severity(severity: int = GENIE_LOG_WARNING) -> None:
    if severity < GENIE_LOG_DEBUG or severity > GENIE_LOG_ERROR:
        raise ValueError("GENIE log severity must be 0 DEBUG, 1 INFO, 2 WARNING, or 3 ERROR")
    _check(genie.GenieSetLogSeverity(severity))


set_log_severity(GENIE_LOG_WARNING)


def _thread_count(threads: int) -> int:
    requested = max(1, int(threads))
    max_threads = max(1, int(os.getenv("DATAPLUG_MGB_MAX_THREADS", "1")))
    return min(requested, max_threads)


def _read_returned_buffer(ptr: C.c_void_p, size: C.c_uint64) -> bytes:
    if not ptr.value or size.value == 0:
        return b""
    remote = cast(ptr.value, C.c_char * int(size.value), cdll=genie)
    return string_at(remote, int(size.value))


def access_unit_count(path: str | Path) -> int:
    count = C.c_uint64()
    _check(genie.GenieGetAccessUnitCount(_path(path), byref(count, cdll=genie)))
    return int(count.value)


def decompress_access_unit(
    path: str | Path,
    access_unit_id: int,
    *,
    reference_file: str | Path | None = None,
    working_dir: str | Path | None = None,
    threads: int = 1,
) -> bytes:
    data = C.c_void_p()
    size = C.c_uint64()
    _check(genie.GenieDecompressAccessUnitToFastq(
        _path(path),
        access_unit_id,
        _optional_path(reference_file),
        _optional_path(working_dir),
        _thread_count(threads),
        byref(data, cdll=genie),
        byref(size, cdll=genie),
    ))
    try:
        return _read_returned_buffer(data, size)
    finally:
        if data.value:
            genie.GenieFree(data.value)


def decompress_access_unit_to_file(
    path: str | Path,
    access_unit_id: int,
    output_file: str | Path,
    *,
    reference_file: str | Path | None = None,
    working_dir: str | Path | None = None,
    threads: int = 1,
) -> None:
    _check(genie.GenieDecompressAccessUnit(
        _path(path),
        access_unit_id,
        _path(output_file),
        _optional_path(reference_file),
        _optional_path(working_dir),
        _thread_count(threads),
    ))


def scan_access_unit(
    path: str | Path,
    access_unit_id: int,
    *,
    reference_file: str | Path | None = None,
    working_dir: str | Path | None = None,
    threads: int = 1,
) -> int:
    if working_dir is None:
        work_path = Path.cwd()
    else:
        work_path = Path(working_dir)
        work_path.mkdir(parents=True, exist_ok=True)
    output_file = work_path / f"access_unit_{access_unit_id}.fastq"
    decompress_access_unit_to_file(
        path,
        access_unit_id,
        output_file,
        reference_file=reference_file,
        working_dir=working_dir,
        threads=threads,
    )
    size = output_file.stat().st_size
    output_file.unlink()
    return size
