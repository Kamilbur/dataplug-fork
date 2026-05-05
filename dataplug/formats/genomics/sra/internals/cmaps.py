from __future__ import annotations

import ctypes
import logging
import os
import struct
import sys
from ctypes import (
    POINTER,
    Structure,
    c_char_p,
    c_int,
    c_size_t,
    c_uint64,
    c_void_p,
    pointer,
    sizeof,
)
from multiprocessing.managers import SharedMemoryManager
from typing import TYPE_CHECKING, TypeAlias

from cdlml import get_var

logger = logging.getLogger(__name__)


libvdb_path = os.environ.get("NCBI_VDB_SO_PATH", "")


_endianness = "<" if sys.byteorder == "little" else ">"
_uint64_t_fmt = _endianness + "Q"
_size_t_fmt = _endianness
if sys.maxsize < 2**32:
    _size_t_fmt += "I"
else:
    _size_t_fmt += "Q"

uint64_t_size = sizeof(c_uint64)


class ShmInfo(Structure):
    _fields_ = (
        ("fd", c_int),
        ("ptr", c_void_p),
        ("length", c_size_t),
    )


if TYPE_CHECKING:
    from ctypes import _Pointer

    ShmInfo_p: TypeAlias = _Pointer[ShmInfo]
else:
    ShmInfo_p = POINTER(ShmInfo)

_ncbi_vdb_mapping = None


class NcbiVdbMapping:
    def __init__(self, libvdb_path: str) -> None:
        self._lib = ctypes.cdll.LoadLibrary(libvdb_path)
        try:
            self._lib.register_shmem.restype = None
            self._lib.register_shmem.argtypes = [ShmInfo_p, c_char_p]

            self._lib.close_shmem.restype = None
            self._lib.close_shmem.argtypes = [ShmInfo_p]

            self._dp_sra_size = get_var(self._lib, c_size_t, "dp_sra_size")
            self._dp_mode = get_var(self._lib, c_int, "dp_mode")

            self.shm_buf = get_var(self._lib, ShmInfo, "shm_buf")
            self.mmap_buf = get_var(self._lib, ShmInfo, "mmap_buf")
            self.pread_buf = get_var(self._lib, ShmInfo, "pread_buf")
        except AttributeError as e:
            msg = f"Probably NCBI_VDB_SO_PATH environment variable is not set.\n[{e!s}]"
            raise AttributeError(msg) from e

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
        self._lib.register_shmem(shm_p, name)

    def close_shmem(self, shm_p: ShmInfo_p) -> None:
        self._lib.close_shmem(shm_p)


def get_ncbi_vdb_mapping(libvdb_path: str | None = None) -> NcbiVdbMapping:
    global _ncbi_vdb_mapping
    if _ncbi_vdb_mapping is None:
        if libvdb_path is None:
            msg = "libvdb_path must be provided on first call"
            raise ValueError(msg)
        _ncbi_vdb_mapping = NcbiVdbMapping(libvdb_path)
    return _ncbi_vdb_mapping


def unpack_from(buf, pos, size=uint64_t_size, fmt=_uint64_t_fmt):
    return struct.unpack(fmt, buf[pos : pos + size])[0]


def pack_to(buf, pos, val, size=uint64_t_size, fmt=_uint64_t_fmt):
    buf[pos : pos + size] = struct.pack(fmt, val)


def ibuf_to_list(buf, length, total_size):
    retval = []

    idx = 0
    for _ in range(length // 2):
        pos = unpack_from(buf, idx)
        size = unpack_from(buf, idx + sizeof(c_uint64), size=sizeof(c_size_t), fmt=_size_t_fmt)
        # Additional 4 bytes for alignment purposes
        if pos + size + 4 < total_size:
            size += 4
        retval.append((pos, pos + size))
        idx += sizeof(c_uint64) + sizeof(c_size_t)

    return retval


def save_mmaps(mmaps, buf, total_size):
    nmmaps = unpack_from(buf, 0)
    if nmmaps == 0:
        return
    tuples = ibuf_to_list(buf[sizeof(c_uint64) :], 2 * nmmaps, total_size)
    mmaps.extend(tuples)
    pack_to(buf, 0, 0)


def save_preads(preads, buf, idx, total_size):
    npreads = unpack_from(buf, 0)
    if npreads == 0:
        return
    tuples = ibuf_to_list(buf[sizeof(c_uint64) :], 2 * npreads, total_size)
    preads[idx] = tuples
    pack_to(buf, 0, 0)


class NcbiVdbSharedMemory:
    def __init__(self, shm_info: ShmInfo, size: int):
        self.size = size
        self.shm_info = shm_info
        self.shared_memory_manager = SharedMemoryManager()
        self.mapping = get_ncbi_vdb_mapping(libvdb_path)

    @property
    def buf(self) -> memoryview[int]:
        buf = self.shared_memory.buf
        if buf is None:
            msg = "Shared memory is not initialized."
            raise RuntimeError(msg)
        return buf

    @property
    def name(self) -> str:
        return self.shared_memory.name

    def __enter__(self):
        self._init_shared_memory()
        self._init_shm_info()
        return self.shared_memory

    def __exit__(self, type, value, traceback):
        self.mapping.close_shmem(pointer(self.shm_info))
        self.shared_memory_manager.shutdown()

    def _init_shared_memory(self):
        self.shared_memory_manager.start()
        self.shared_memory = self.shared_memory_manager.SharedMemory(size=self.size)

    def _init_shm_info(self):
        self.shm_info.length = self.size
        self.mapping.register_shmem(pointer(self.shm_info), self.name.encode())
