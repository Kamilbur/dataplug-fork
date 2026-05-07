from __future__ import annotations

import ctypes as C

import cdlml


def _func(name: str):
    return getattr(cdlml, name, None) or getattr(C, name)


def _is_carg(obj) -> bool:
    return type(obj).__name__ == "CArgObject"


def _contextual(name: str, obj, *extra, cdll=None):
    fn = getattr(cdlml, name, None)
    if fn is not None:
        attempts = [(obj, *extra)]
        if cdll is not None:
            attempts.extend(((cdll, obj, *extra), (obj, *extra, cdll)))
        fallback = None
        for args in attempts:
            try:
                result = fn(*args)
            except TypeError:
                continue
            fallback = result
            if cdll is None or hasattr(result, "_addr") or not _is_carg(result):
                return result
        if fallback is not None:
            return fallback
    return getattr(C, name)(obj, *extra)


def _addr(obj):
    try:
        return addressof(obj)
    except TypeError:
        return obj


def addressof(obj):
    return _func("addressof")(obj)


def byref(obj, cdll=None):
    fn = getattr(cdlml, "byref", None)
    if fn is not None and cdll is not None:
        try:
            result = fn(obj, cdll=cdll)
        except TypeError:
            pass
        else:
            if hasattr(result, "_addr") or not _is_carg(result):
                return result
    return _contextual("byref", obj, cdll=cdll)


def cast(obj, typ, cdll=None):
    fn = getattr(cdlml, "cast", None)
    if fn is not None and cdll is not None:
        try:
            result = fn(obj, typ, cdll=cdll)
        except TypeError:
            pass
        else:
            return result
    try:
        return _contextual("cast", obj, typ, cdll=cdll)
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


def pointer_for(obj, cdll=None):
    fn = getattr(cdlml, "pointer", None)
    if fn is not None and cdll is not None:
        try:
            return fn(obj, cdll=cdll)
        except TypeError:
            pass
    return pointer(obj)


def string_at(ptr, size=-1):
    return _func("string_at")(ptr, size)


def sizeof(obj):
    return _func("sizeof")(obj)
