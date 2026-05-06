from __future__ import annotations

import os
import shutil
import subprocess
import tarfile
import urllib.request
from pathlib import Path

from setuptools import setup
from setuptools.command.build_py import build_py as _build_py


GENIE_REPO = "https://github.com/Kamilbur/genie.git"
GENIE_REF = os.environ.get("MGB_NATIVE_GENIE_REF", os.environ.get("GENIE_NATIVE_GENIE_REF", ""))

ZLIB_REPO = "https://github.com/madler/zlib.git"
ZLIB_REF = "v1.3.1"

BZIP2_REPO = "https://sourceware.org/git/bzip2.git"
BZIP2_REF = "bzip2-1.0.8"

XZ_REF = "5.4.6"
XZ_ARCHIVE = f"https://github.com/tukaani-project/xz/releases/download/v{XZ_REF}/xz-{XZ_REF}.tar.gz"

ZSTD_REPO = "https://github.com/facebook/zstd.git"
ZSTD_REF = "v1.5.6"

LIBBSC_REPO = "https://github.com/IlyaGrebnov/libbsc.git"
LIBBSC_REF = "v3.3.12"

HTSLIB_REF = "1.20"
HTSLIB_ARCHIVE = (
    f"https://downloads.sourceforge.net/project/samtools/samtools/{HTSLIB_REF}/htslib-{HTSLIB_REF}.tar.bz2"
)


_HERE = Path(__file__).parent.resolve()
_ROOT = _HERE.parents[1]
_WORK = _HERE / "_build"
_DEPS = _WORK / "dependencies"
_INSTALL = _WORK / "install"
_INTERNALS = _ROOT / "dataplug" / "formats" / "genomics" / "mgb" / "internals"


class BuildPy(_build_py):
    def run(self):
        _build_dependencies()
        _build_genie()
        _build_shims()
        super().run()


def _jobs() -> str:
    return str(max(1, (os.cpu_count() or 2) // 2))


def _run(args, cwd: Path | None = None, env: dict[str, str] | None = None):
    print("mgb_native:", " ".join(map(str, args)))
    subprocess.run([str(arg) for arg in args], cwd=cwd, env=env, check=True)


def _clone(repo: str, ref: str, dst: Path, *, submodules: bool = False):
    if dst.exists():
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["git", "clone", "--depth", "1"]
    if ref:
        cmd.extend(["--branch", ref])
    cmd.extend([repo, dst])
    _run(cmd)
    if submodules:
        _run(["git", "submodule", "update", "--init", "--recursive"], cwd=dst)


def _download_archive(url: str, dst: Path):
    if dst.exists():
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    print(f"mgb_native: downloading {url}")
    urllib.request.urlretrieve(url, dst)


def _extract_tarball(archive: Path, dst: Path):
    if dst.exists():
        return
    tmp = dst.parent / f".{dst.name}.tmp"
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir(parents=True)
    with tarfile.open(archive, "r:*") as tf:
        tf.extractall(tmp)
    matches = [path.parent for path in tmp.rglob("configure") if path.is_file()]
    if not matches:
        raise RuntimeError(f"configure not found in {archive}")
    extracted = matches[0]
    shutil.move(str(extracted), dst)
    shutil.rmtree(tmp)


def _env() -> dict[str, str]:
    env = os.environ.copy()
    include = _INSTALL / "include"
    lib_dirs = _lib_dirs()
    pkg_config_dirs = [lib / "pkgconfig" for lib in lib_dirs]
    env["CPPFLAGS"] = f"-I{include} " + env.get("CPPFLAGS", "")
    env["CFLAGS"] = f"-fPIC -I{include} " + env.get("CFLAGS", "")
    env["CXXFLAGS"] = f"-fPIC -I{include} -I{include / 'libbsc'} " + env.get("CXXFLAGS", "")
    env["LDFLAGS"] = " ".join(f"-L{lib} -Wl,-rpath,{lib}" for lib in lib_dirs) + " " + env.get("LDFLAGS", "")
    env["PKG_CONFIG_PATH"] = ":".join(str(path) for path in pkg_config_dirs) + f":{env.get('PKG_CONFIG_PATH', '')}"
    env["LD_LIBRARY_PATH"] = ":".join(str(path) for path in lib_dirs) + f":{env.get('LD_LIBRARY_PATH', '')}"
    return env


def _primary_lib_dir() -> Path:
    return _INSTALL / "lib"


def _lib_dirs() -> list[Path]:
    return [_INSTALL / "lib", _INSTALL / "lib64"]


def _lib(name: str) -> Path:
    for lib_dir in _lib_dirs():
        for candidate in lib_dir.glob(name):
            return candidate
    raise RuntimeError(f"required library not found: {name} in {_lib_dirs()}")


def _lib_exists(name: str) -> bool:
    return any(lib_dir.joinpath(name).exists() for lib_dir in _lib_dirs())


def _copy_first_lib(pattern: str, dst: Path, missing_msg: str) -> None:
    for lib_dir in _lib_dirs():
        libs = sorted(lib_dir.glob(pattern))
        if libs:
            shutil.copy2(libs[-1], dst)
            return
    raise RuntimeError(missing_msg)


def _plain_so_name(path: Path) -> str:
    name = path.name
    marker = ".so"
    idx = name.find(marker)
    if idx == -1:
        return name
    return name[: idx + len(marker)]


def _build_dependencies():
    _DEPS.mkdir(parents=True, exist_ok=True)
    _INSTALL.mkdir(parents=True, exist_ok=True)

    _build_zlib()
    _build_bzip2()
    _build_xz()
    _build_zstd()
    _build_libbsc()
    _build_htslib()


def _build_zlib():
    src = _DEPS / "zlib"
    build = _WORK / "build-zlib"
    if _lib_exists("libz.so"):
        return
    _clone(ZLIB_REPO, ZLIB_REF, src)
    _run([
        "cmake",
        "-S", src,
        "-B", build,
        "-DCMAKE_BUILD_TYPE=Release",
        f"-DCMAKE_INSTALL_PREFIX={_INSTALL}",
        "-DCMAKE_INSTALL_LIBDIR=lib",
        "-DCMAKE_POSITION_INDEPENDENT_CODE=ON",
        "-DBUILD_SHARED_LIBS=ON",
    ], env=_env())
    _run(["cmake", "--build", build, "--parallel", _jobs()], env=_env())
    _run(["cmake", "--install", build], env=_env())


def _build_bzip2():
    src = _DEPS / "bzip2"
    marker = _primary_lib_dir() / "libbz2.so"
    if _lib_exists("libbz2.so"):
        return
    _clone(BZIP2_REPO, BZIP2_REF, src)
    _run(["make", "-f", "Makefile-libbz2_so", "-j", _jobs()], cwd=src, env=_env())
    (_INSTALL / "include").mkdir(parents=True, exist_ok=True)
    _primary_lib_dir().mkdir(parents=True, exist_ok=True)
    shutil.copy2(src / "bzlib.h", _INSTALL / "include" / "bzlib.h")
    for so in src.glob("libbz2.so*"):
        shutil.copy2(so, _primary_lib_dir() / so.name)
    if not marker.exists():
        _copy_first_lib("libbz2.so.*", marker, "libbz2.so not built")


def _build_xz():
    src = _DEPS / "xz"
    if _lib_exists("liblzma.so"):
        return
    if src.exists() and not (src / "configure").exists():
        shutil.rmtree(src)
    archive = _DEPS / f"xz-{XZ_REF}.tar.gz"
    _download_archive(XZ_ARCHIVE, archive)
    _extract_tarball(archive, src)
    _run([
        "./configure",
        f"--prefix={_INSTALL}",
        "--disable-static",
        "--enable-shared",
        "--disable-doc",
        "--disable-nls",
    ], cwd=src, env=_env())
    _run(["make", "-j", _jobs()], cwd=src, env=_env())
    _run(["make", "install"], cwd=src, env=_env())


def _build_zstd():
    src = _DEPS / "zstd"
    if _lib_exists("libzstd.so"):
        return
    _clone(ZSTD_REPO, ZSTD_REF, src)
    _run([
        "make",
        "-j", _jobs(),
        f"PREFIX={_INSTALL}",
        "BUILD_SHARED=1",
        "BUILD_STATIC=0",
        "install",
    ], cwd=src, env=_env())


def _build_libbsc():
    src = _DEPS / "libbsc"
    build = _WORK / "build-libbsc"
    marker = _primary_lib_dir() / "libbsc.so"
    if _lib_exists("libbsc.so"):
        return
    _clone(LIBBSC_REPO, LIBBSC_REF, src)
    _run([
        "cmake",
        "-S", src,
        "-B", build,
        "-DCMAKE_BUILD_TYPE=Release",
        f"-DCMAKE_INSTALL_PREFIX={_INSTALL}",
        "-DCMAKE_INSTALL_LIBDIR=lib",
        "-DCMAKE_POSITION_INDEPENDENT_CODE=ON",
        "-DBUILD_SHARED_LIBS=ON",
        "-DBSC_BUILD_SHARED_LIB=ON",
    ], env=_env())
    _run(["cmake", "--build", build, "--parallel", _jobs()], env=_env())
    _run(["cmake", "--install", build], env=_env())

    include_dir = _INSTALL / "include" / "libbsc"
    include_dir.mkdir(parents=True, exist_ok=True)
    header = next(_INSTALL.rglob("libbsc.h"), None) or next(src.rglob("libbsc.h"))
    header_dst = include_dir / "libbsc.h"
    if header.resolve() != header_dst.resolve():
        shutil.copy2(header, header_dst)
    if not marker.exists():
        _copy_first_lib("libbsc.so*", marker, "libbsc.so not built")


def _build_htslib():
    src = _DEPS / "htslib"
    if _lib_exists("libhts.so"):
        return
    if src.exists() and not (src / "configure").exists():
        shutil.rmtree(src)
    archive = _DEPS / f"htslib-{HTSLIB_REF}.source.tar.bz2"
    _download_archive(HTSLIB_ARCHIVE, archive)
    _extract_tarball(archive, src)
    _run([
        "./configure",
        f"--prefix={_INSTALL}",
        "--disable-libcurl",
        "--disable-s3",
        "--disable-gcs",
        "--without-libdeflate",
        "--enable-shared",
        "--disable-static",
    ], cwd=src, env=_env())
    _run(["make", "-j", _jobs()], cwd=src, env=_env())
    _run(["make", "install"], cwd=src, env=_env())


def _build_genie():
    src = _WORK / "genie"
    build = _WORK / "genie-build"
    _clone(GENIE_REPO, GENIE_REF, src)
    _patch_genie_source(src)
    build.mkdir(parents=True, exist_ok=True)

    install_include = _INSTALL / "include"
    install_lib = _primary_lib_dir()
    lib_rpath = ";".join(str(path) for path in _lib_dirs())
    cmake_args = [
        "cmake",
        "-S", src,
        "-B", build,
        "-DCMAKE_BUILD_TYPE=Release",
        f"-DCMAKE_PREFIX_PATH={_INSTALL}",
        f"-DCMAKE_BUILD_RPATH=$ORIGIN;{lib_rpath}",
        f"-DCMAKE_INSTALL_RPATH=$ORIGIN",
        "-DGENIE_USE_OPENMP=OFF",
        "-DGENIE_SAM_SUPPORT=ON",
        "-DBUILD_TESTS=OFF",
        f"-DBSC_LIBRARY={_lib('libbsc.so')}",
        f"-DBSC_INCLUDE_DIR={install_include / 'libbsc'}",
        f"-DZSTD_LIBRARY={_lib('libzstd.so')}",
        f"-DZSTD_INCLUDE_DIR={install_include}",
        f"-DLZMA_LIBRARY={_lib('liblzma.so')}",
        # This GENIE fork's FindLZMA.cmake probes libbsc.h, not lzma.h.
        f"-DLZMA_INCLUDE_DIR={install_include / 'libbsc'}",
        f"-DHTSLIB_INCLUDE_DIR={install_include}",
        f"-DHTSlib_LIBRARY={_lib('libhts.so')}",
        f"-DZLIB_INCLUDE_DIR={install_include}",
        f"-DZLIB_LIBRARY={_lib('libz.so')}",
        f"-DZLIB_LIBRARIES={_lib('libz.so')}",
        f"-DBZIP2_INCLUDE_DIR={install_include}",
        f"-DBZIP2_LIBRARIES={_lib('libbz2.so')}",
        f"-DCMAKE_CXX_FLAGS=-I{install_include} -I{install_include / 'libbsc'}",
    ]
    _run(cmake_args, env=_env())
    _run(["cmake", "--build", build, "--target", "genie-shared", "--parallel", _jobs()], env=_env())

    libgenie = build / "lib" / "libgenie.so"
    if not libgenie.exists():
        raise RuntimeError("libgenie.so not found after GENIE build")

    _INTERNALS.mkdir(parents=True, exist_ok=True)
    shutil.copy2(libgenie, _INTERNALS / "libgenie.so")
    for pattern in ("libbsc.so*", "libzstd.so*", "liblzma.so*", "libhts.so*", "libz.so*", "libbz2.so*"):
        for lib_dir in _lib_dirs():
            for lib in lib_dir.glob(pattern):
                if lib.is_file():
                    shutil.copy2(lib, _INTERNALS / _plain_so_name(lib))
    print(f"mgb_native: installed libgenie.so and dependency libs -> {_INTERNALS}")


def _patch_genie_source(src: Path):
    src_cmake = src / "src" / "CMakeLists.txt"
    lines = src_cmake.read_text(encoding="utf-8").splitlines()
    patched = "\n".join(line for line in lines if line.strip() != "add_subdirectory(apps)") + "\n"
    text = "\n".join(lines) + "\n"
    if patched != text:
        src_cmake.write_text(patched, encoding="utf-8")


def _build_shims():
    shims_c = _INTERNALS / "shims.c"
    if not shims_c.exists():
        raise RuntimeError(f"shims.c not found at {shims_c}")
    _run([
        "gcc",
        "-shared",
        "-fPIC",
        "-o", _INTERNALS / "libmgbshims.so",
        shims_c,
        "-ldl",
        "-lpthread",
        "-Wall",
        "-Wextra",
        "-O2",
    ])
    print(f"mgb_native: built libmgbshims.so -> {_INTERNALS}")


setup(
    name="dataplug-mgb-native",
    version="1.0.0",
    packages=["mgb_native"],
    cmdclass={"build_py": BuildPy},
)
