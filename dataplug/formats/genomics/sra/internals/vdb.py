import ctypes as C
import importlib.machinery
import importlib.util
import sys
import sysconfig
import types
from pathlib import Path


def _install_cdlml_extension_shim() -> None:
    """Allow cdlml to use dlmopen when its package shim is stale.

    Some editable cdlml installs have a working `_cdlml._cdlml` extension but
    an `_cdlml/__init__.py` that imports `_is_available`, while the built
    extension does not export it. cdlml treats that ImportError as "no glibc
    dlmopen" and falls back to the subprocess server, which is not equivalent
    for VDB interposition. Install a tiny module shim exposing the symbols
    cdlml needs when the extension is otherwise loadable.
    """
    try:
        import _cdlml  # noqa: F401

        return
    except ImportError:
        pass

    suffix = sysconfig.get_config_var("EXT_SUFFIX") or ".so"
    for base in map(Path, sys.path):
        candidate = base / "_cdlml" / f"_cdlml{suffix}"
        if not candidate.exists():
            continue

        pkg = types.ModuleType("_cdlml")
        pkg.__path__ = [str(candidate.parent)]
        sys.modules["_cdlml"] = pkg

        loader = importlib.machinery.ExtensionFileLoader("_cdlml._cdlml", str(candidate))
        spec = importlib.util.spec_from_file_location("_cdlml._cdlml", candidate, loader=loader)
        if spec is None:
            continue
        mod = importlib.util.module_from_spec(spec)
        sys.modules["_cdlml._cdlml"] = mod
        try:
            loader.exec_module(mod)
        except ImportError:
            sys.modules.pop("_cdlml._cdlml", None)
            sys.modules.pop("_cdlml", None)
            continue

        pkg._dlmopen = mod._dlmopen
        pkg._dlmstop = mod._dlmstop
        pkg._is_available = getattr(mod, "_is_available", lambda: True)
        return


_install_cdlml_extension_shim()

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
