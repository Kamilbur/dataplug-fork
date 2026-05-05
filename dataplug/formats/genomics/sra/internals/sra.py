import contextlib
import ctypes as C

from .vdb import vdb
from .cdlml_compat import byref
from .vdb_types import to_char_p

SRALines = list[tuple[str, ...]]


_COLUMN_NAMES = (
    "READ",
    "(INSDC:quality:text:phred_33)QUALITY",
    "NAME",
    "READ_START",
    "READ_LEN",
)


def _release(name, handle):
    if handle is not None and handle.value:
        with contextlib.suppress(Exception):
            getattr(vdb, name)(handle)


def _open_table(mgr, path):
    tab = C.c_void_p()
    if vdb.VDBManagerOpenTableRead(mgr, byref(tab), C.c_void_p(0), to_char_p(path)) == 0:
        return tab, None
    db = C.c_void_p()
    if vdb.VDBManagerOpenDBRead(mgr, byref(db), C.c_void_p(0), to_char_p(path)) != 0:
        raise ValueError("Not an SRA-object: " + path)
    if vdb.VDatabaseOpenTableRead(db, byref(tab), to_char_p("SEQUENCE")) == 0:
        return tab, db
    _release("VDatabaseRelease", db)
    raise ValueError("Not an SRA-object: " + path)


def _open_cursor(path):
    nat_dir = C.c_void_p()
    mgr = C.c_void_p()
    cur = C.c_void_p()
    vdb.KDirectoryNativeDir_v1(byref(nat_dir))
    vdb.VDBManagerMakeRead(byref(mgr), nat_dir)
    tab, db = _open_table(mgr, path)
    vdb.VTableCreateCursorRead(tab, byref(cur))
    return cur, tab, db, mgr, nat_dir


class VColumns:
    def __init__(self, handles):
        from .vdb_types import VColumn

        cur, tab, db, mgr, nat_dir = handles
        self.cur = cur
        self._tab = tab
        self._db = db
        self._mgr = mgr
        self._nat_dir = nat_dir
        self._closed = False
        self.columns = []
        for name in _COLUMN_NAMES:
            idx = C.c_int()
            vdb.VCursorAddColumn(cur, byref(idx), to_char_p(name))
            self.columns.append(VColumn(idx.value, cur))
        vdb.VCursorOpen(cur)
        for col in self.columns:
            col.update()

    def __len__(self):
        return self.columns[0].row_range()[1]

    @classmethod
    def from_filepath(cls, path):
        return cls(_open_cursor(path))

    def close(self):
        if self._closed:
            return
        self._closed = True
        _release("VCursorRelease", self.cur)
        _release("VTableRelease", self._tab)
        _release("VDatabaseRelease", self._db)
        _release("VDBManagerRelease", self._mgr)
        _release("KDirectoryRelease_v1", self._nat_dir)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()

    def __del__(self):
        self.close()


def _format_spot(acc, idx, read_str, qual_str, name, starts, lengths, split=True):
    header = f"{acc}.{idx} {name} length={len(read_str)}"
    if not split:
        return (f"@{header}\n{read_str}\n+{header}\n{qual_str}\n",)
    result = []
    for start, length in zip(starts, lengths, strict=True):
        if length > 0:
            r = read_str[start : start + length]
            q = qual_str[start : start + length]
            result.append(f"@{header}\n{r}\n+{header}\n{q}\n")
    return tuple(result)
