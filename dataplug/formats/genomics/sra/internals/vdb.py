from pathlib import Path
import ctypes as C

from cdlml import Preload, PreloadedCDLL

_HERE = Path(__file__).parent

preloads = [
    Preload(name=b"libshims.so", path=str(_HERE / "libshims.so").encode(), alias="shims"),
]

vdb = PreloadedCDLL(str(_HERE / "libncbi-vdb.so"), preloads=preloads)


def _prototype(name, restype, argtypes):
    func = getattr(vdb, name)
    func.restype = restype
    func.argtypes = argtypes


_P = C.c_void_p
_prototype("KDirectoryNativeDir_v1", C.c_int, [_P])
_prototype("VDBManagerMakeRead", C.c_int, [_P, _P])
_prototype("VDBManagerOpenTableRead", C.c_int, [_P, _P, _P, C.c_char_p])
_prototype("VDBManagerOpenDBRead", C.c_int, [_P, _P, _P, C.c_char_p])
_prototype("VDatabaseOpenTableRead", C.c_int, [_P, _P, C.c_char_p])
_prototype("VTableCreateCursorRead", C.c_int, [_P, _P])
_prototype("VCursorAddColumn", C.c_int, [_P, _P, C.c_char_p])
_prototype("VCursorOpen", C.c_int, [_P])
_prototype("VCursorDatatype", C.c_int, [_P, C.c_int, _P, _P])
_prototype("VCursorCellDataDirect", C.c_int, [_P, C.c_longlong, C.c_int, _P, _P, _P, _P])
_prototype("VCursorIdRange", C.c_int, [_P, C.c_int, _P, _P])
for _name in (
    "VCursorRelease",
    "VTableRelease",
    "VDatabaseRelease",
    "VDBManagerRelease",
    "KDirectoryRelease_v1",
):
    _prototype(_name, C.c_int, [_P])
