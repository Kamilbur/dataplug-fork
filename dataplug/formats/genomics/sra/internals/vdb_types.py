import ctypes as C
from dataclasses import dataclass, field
from enum import Enum

from .cdlml_compat import byref, cast, string_at
from . import vdb as vdb_module


def to_char_p(s):
    if isinstance(s, str):
        return s.encode()
    if isinstance(s, C.c_char_p):
        return s.value
    return s


class Typedecl(C.Structure):
    _fields_ = (
        ("type_id", C.c_int),
        ("dim", C.c_int),
    )


class Typedesc(C.Structure):
    _fields_ = (
        ("bits", C.c_int),
        ("dim", C.c_int),
        ("domain", C.c_int),
    )


uint_xf = {1: C.c_ubyte, 8: C.c_ubyte, 16: C.c_ushort, 32: C.c_uint, 64: C.c_ulonglong}

int_xf = {1: C.c_byte, 8: C.c_byte, 16: C.c_short, 32: C.c_int, 64: C.c_longlong}

float_xf = {32: C.c_float, 64: C.c_double}

txt_xf = {8: C.c_char}


class TypeDomain(Enum):
    Bool = 1
    UInt = 2
    Int = 3
    Float = 4
    Ascii = 5
    Unicode = 6


type_xf = {
    TypeDomain.Int.value: (int_xf, C.c_byte),
    TypeDomain.Bool.value: (uint_xf, C.c_ubyte),
    TypeDomain.UInt.value: (uint_xf, C.c_ubyte),
    TypeDomain.Float.value: (float_xf, C.c_double),
    TypeDomain.Ascii.value: (txt_xf, C.c_char),
    TypeDomain.Unicode.value: (txt_xf, C.c_char),
}


@dataclass
class VColumn:
    idx: int
    cur: C.c_void_p
    cast: str = ""
    name: str = ""
    column_type: tuple | None = None
    tdec: Typedecl = field(default_factory=lambda: Typedecl(0, 0))
    tdes: Typedesc = field(default_factory=lambda: Typedesc(0, 0, 0))

    def cast_name(self, name: str) -> str:
        p1 = name.find("(")
        p2 = name.find(")")
        if p1 > -1 and p2 > -1:
            self.cast = name[p1 + 1 : p2]
            self.name = name[p2 + 1 :]
        else:
            self.name = name

    def update(self):
        vdb = vdb_module.vdb
        vdb.VCursorDatatype(self.cur, self.idx, byref(self.tdec, vdb), byref(self.tdes, vdb))
        (dict, dflt) = type_xf[self.tdes.domain]
        self.column_type = dict[self.tdes.bits]
        if self.column_type is None:
            self.column_type = dflt

    def read(self, row):
        vdb = vdb_module.vdb
        row_id = C.c_longlong(row)
        elem_bits = C.c_int()
        data = C.c_void_p()
        row_len = C.c_int()
        vdb.VCursorCellDataDirect(
            self.cur,
            row_id.value,
            self.idx,
            byref(elem_bits, vdb),
            byref(data, vdb),
            0,
            byref(row_len, vdb),
        )
        if self.column_type == C.c_char:
            tmp = string_at(cast(data.value or 0, C.POINTER(C.c_char), vdb), row_len.value)
            if isinstance(tmp, bytes):
                return tmp.decode("utf-8")
            return tmp
        e_count = row_len.value
        if elem_bits.value < 8:
            e_count *= elem_bits.value
            e_count //= 8
        typed_ptr = cast(data.value or 0, C.POINTER(self.column_type), vdb)
        raw = string_at(typed_ptr, e_count * C.sizeof(self.column_type))
        arr_type = self.column_type * e_count
        return list(arr_type.from_buffer_copy(raw))

    def row_range(self):
        vdb = vdb_module.vdb
        first = C.c_longlong()
        count = C.c_longlong()
        vdb.VCursorIdRange(self.cur, self.idx, byref(first, vdb), byref(count, vdb))
        return (first.value, count.value)
