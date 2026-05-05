import subprocess
import os
import shutil
from pathlib import Path
from setuptools import setup
from setuptools.command.build_py import build_py as _build_py

NCBI_VDB_REPO = "https://github.com/ncbi/ncbi-vdb"
#NCBI_VDB_REPO = "https://github.com/Kamilbur/ncbi-vdb-dp/tree/optimized#"
#NCBI_VDB_REPO = "https://github.com/Kamilbur/ncbi-vdb-dp"

_HERE = Path(__file__).parent
DEST = _HERE.parent / "dataplug" / "formats" / "genomics" / "sra"


class BuildPy(_build_py):
    def run(self):
        _build_ncbi_vdb()
        super().run()


def _build_ncbi_vdb():
    work = _HERE / "_build"
    repo = work / "ncbi-vdb"
    cmake_build = work / "cmake-build"

    work.mkdir(exist_ok=True)
    if not repo.exists():
        subprocess.run(
            # ["git", "clone", "-b", "optimized", "--depth", "1", NCBI_VDB_REPO, str(repo)],
            ["git", "clone", "--depth", "1", NCBI_VDB_REPO, str(repo)],
            check=True,
        )

    cmake_build.mkdir(exist_ok=True)
    subprocess.run(["cmake", "-S", str(repo), "-B", str(cmake_build)], check=True)
    subprocess.run(["cmake", "--build", str(cmake_build), "--parallel", str(max(1, os.cpu_count() // 2))], check=True)

    so_files = list(cmake_build.rglob("libncbi-vdb.so"))
    if not so_files:
        raise RuntimeError("libncbi-vdb.so not found after build")

    DEST.mkdir(parents=True, exist_ok=True)
    shutil.copy2(so_files[0], DEST / "libncbi-vdb.so")
    shutil.copy2(repo / "py_vdb" / "vdb.py", DEST / "vdb.py")
    print(f"sra_native: installed libncbi-vdb.so and vdb.py → {DEST}")


setup(
    name="dataplug-sra-native",
    version="1.0.0",
    packages=["sra_native"],
    cmdclass={"build_py": BuildPy},
)
