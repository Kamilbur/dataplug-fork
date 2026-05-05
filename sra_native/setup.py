import os
import shutil
import subprocess
from pathlib import Path

from setuptools import setup
from setuptools.command.build_py import build_py as _build_py

NCBI_VDB_REPO = "https://github.com/ncbi/ncbi-vdb"

_HERE = Path(__file__).parent
_INTERNALS = _HERE.parent / "dataplug" / "formats" / "genomics" / "sra" / "internals"


class BuildPy(_build_py):
    def run(self):
        _build_ncbi_vdb()
        _build_shims()
        super().run()


def _build_ncbi_vdb():
    work = _HERE / "_build"
    repo = work / "ncbi-vdb"
    cmake_build = work / "cmake-build"

    work.mkdir(exist_ok=True)
    if not repo.exists():
        subprocess.run(
            ["git", "clone", "--depth", "1", NCBI_VDB_REPO, str(repo)],
            check=True,
        )

    cmake_build.mkdir(exist_ok=True)
    subprocess.run(
        ["cmake", "-S", str(repo), "-B", str(cmake_build), "-DCMAKE_BUILD_TYPE=Debug"], check=True
    )
    subprocess.run(
        ["cmake", "--build", str(cmake_build), "--parallel", str(max(1, os.cpu_count() // 2))],
        check=True,
    )

    so_files = list(cmake_build.rglob("libncbi-vdb.so"))
    if not so_files:
        raise RuntimeError("libncbi-vdb.so not found after build")

    _INTERNALS.mkdir(parents=True, exist_ok=True)
    shutil.copy2(so_files[0], _INTERNALS / "libncbi-vdb.so")
    print(f"sra_native: installed libncbi-vdb.so → {_INTERNALS}")


def _build_shims():
    shims_c = _INTERNALS / "shims.c"
    if not shims_c.exists():
        raise RuntimeError(f"shims.c not found at {shims_c}")

    subprocess.run(
        [
            "gcc",
            "-shared",
            "-fPIC",
            "-o",
            str(_INTERNALS / "libshims.so"),
            str(shims_c),
            "-ldl",
            "-lpthread",
            "-Wall",
            "-Wextra",
            "-O2",
        ],
        check=True,
    )
    print(f"sra_native: built libshims.so → {_INTERNALS}")


setup(
    name="dataplug-sra-native",
    version="1.0.0",
    packages=["sra_native"],
    cmdclass={"build_py": BuildPy},
)
