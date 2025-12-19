import importlib
import os
import sys

libvdb_path = os.environ.get('NCBI_VDB_SO_PATH', '')
pyvdb_dir_path = 'PY_VDB_DIR'


def _import_lib(libname, libpath=''):
    try:
        return importlib.import_module(libname)
    except ImportError:
        pass
    d = os.getenv(libpath)
    if d and d not in sys.path:
        sys.path.insert(0, d)
        return importlib.import_module(libname)
    msg = f"Could not import {libname} module."
    raise ImportError(msg)


try:
    vdb = _import_lib('vdb', pyvdb_dir_path)
except ImportError as e:
    msg = (
        "Could not import vdb module. "
        "Please make sure NCBI VDB SDK is installed and "
        "PY_VDB_DIR environment variable is set correctly."
    )
    raise ImportError(msg) from e

_vdb_manager = vdb.manager(mode=vdb.OpenMode.Read, path=libvdb_path)


def get_vdb():
    return vdb


def get_vdb_manager():
    return _vdb_manager


def get_table(accession: str):
    try:
        return _vdb_manager.OpenDB(accession).OpenTable("SEQUENCE")
    except vdb.vdb_error as _:
        pass
    try:
        return _vdb_manager.OpenTable(accession)
    except vdb.vdb_error as _:
        pass
    msg = (
        "Could not open sequence table. "
        f"{accession} is not an SRA-object"
    )
    raise ValueError(msg)
