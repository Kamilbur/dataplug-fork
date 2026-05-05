from __future__ import annotations

import ctypes as C

import cdlml


def _func(name: str):
    return getattr(cdlml, name, None) or getattr(C, name)


def _addr(obj):
    try:
        return addressof(obj)
    except TypeError:
        return obj


def addressof(obj):
    return _func("addressof")(obj)


def byref(obj):
    return _func("byref")(obj)


def cast(obj, typ):
    fn = _func("cast")
    try:
        return fn(obj, typ)
    except TypeError:
        return C.cast(obj, typ)


def memmove(dst, src, size):
    fn = _func("memmove")
    try:
        return fn(dst, src, size)
    except TypeError:
        return fn(_addr(dst), _addr(src), size)


def memset(dst, value, size):
    fn = _func("memset")
    try:
        return fn(dst, value, size)
    except TypeError:
        return fn(_addr(dst), value, size)


def pointer(obj):
    return _func("pointer")(obj)


def sizeof(obj):
    return _func("sizeof")(obj)
