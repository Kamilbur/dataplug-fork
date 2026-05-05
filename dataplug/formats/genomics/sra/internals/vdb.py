from pathlib import Path

from cdlml import Preload, PreloadedCDLL

_HERE = Path(__file__).parent

preloads = [
    Preload(name=b"libshims.so", path=str(_HERE / "libshims.so").encode(), alias="shims"),
]

vdb = PreloadedCDLL(str(_HERE / "libncbi-vdb.so"), preloads=preloads)
