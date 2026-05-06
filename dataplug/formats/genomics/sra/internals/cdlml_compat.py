from __future__ import annotations

import ctypes as C

import cdlml


def _func(name: str):
    return getattr(cdlml, name, None) or getattr(C, name)


def _is_carg(obj) -> bool:
    return type(obj).__name__ == "CArgObject"


def _remote_cdll(cdll):
    return getattr(cdll, "_cdll", cdll)


def _contextual(name: str, obj, *extra, cdll=None):
    fn = getattr(cdlml, name, None)
    if fn is not None:
        attempts = [(obj, *extra)]
        if cdll is not None:
            remote = _remote_cdll(cdll)
            attempts.extend(((remote, obj, *extra), (obj, *extra, remote)))
            for keyword in ("cdll", "lib", "library", "server"):
                attempts.append((obj, *extra, {keyword: remote}))
        fallback = None
        for args in attempts:
            try:
                if args and isinstance(args[-1], dict):
                    result = fn(*args[:-1], **args[-1])
                else:
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
    return _contextual("byref", obj, cdll=cdll)


def cast(obj, typ, cdll=None):
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


def sizeof(obj):
    return _func("sizeof")(obj)
