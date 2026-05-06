"""Download the SRR32296234 SRA fixture used by speed benchmarks.

Examples:
    python3 tests/download_srr32296234.py
    python3 tests/download_srr32296234.py --dest /tmp/SRR32296234
    python3 tests/download_srr32296234.py --force
"""

from __future__ import annotations

import argparse
import time
import urllib.error
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DEST = REPO_ROOT / "tests" / "fixtures" / "sra" / "SRR32296234"
URL = "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR32296234/SRR32296234"
EXPECTED_SIZE = 3_279_135_186
CHUNK_SIZE = 8 * 1024 * 1024


def fmt_size(size: int) -> str:
    units = ("B", "KiB", "MiB", "GiB")
    value = float(size)
    for unit in units:
        if value < 1024 or unit == units[-1]:
            return f"{value:.1f} {unit}"
        value /= 1024
    return f"{size} B"


def remote_size(url: str) -> int:
    req = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(req, timeout=30) as resp:
        return int(resp.headers["Content-Length"])


def download(url: str, dest: Path, expected_size: int, force: bool) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")

    if dest.exists() and not force:
        size = dest.stat().st_size
        if size == expected_size:
            print(f"exists ok: {dest} ({fmt_size(size)})")
            return
        raise SystemExit(
            f"{dest} exists with size {size}, expected {expected_size}. "
            "Use --force to replace."
        )

    if force:
        part.unlink(missing_ok=True)

    offset = part.stat().st_size if part.exists() else 0
    if offset > expected_size:
        raise SystemExit(f"partial file too large: {part} ({offset} bytes)")

    headers = {}
    mode = "ab"
    if offset:
        headers["Range"] = f"bytes={offset}-"
        print(f"resuming at {fmt_size(offset)}")
    else:
        print(f"downloading {url}")

    req = urllib.request.Request(url, headers=headers)
    started = time.perf_counter()
    last_report = started
    downloaded = offset

    try:
        with urllib.request.urlopen(req, timeout=60) as resp, part.open(mode) as f:
            status = getattr(resp, "status", None)
            if offset and status != 206:
                raise RuntimeError(f"server did not honor Range request, status={status}")

            while True:
                chunk = resp.read(CHUNK_SIZE)
                if not chunk:
                    break
                f.write(chunk)
                downloaded += len(chunk)

                now = time.perf_counter()
                if now - last_report >= 5:
                    elapsed = max(now - started, 1e-9)
                    speed = (downloaded - offset) / elapsed
                    pct = downloaded / expected_size * 100
                    print(
                        f"{pct:6.2f}% {fmt_size(downloaded)} / "
                        f"{fmt_size(expected_size)} at {fmt_size(int(speed))}/s",
                        flush=True,
                    )
                    last_report = now
    except urllib.error.URLError as exc:
        raise SystemExit(f"download failed: {exc}") from exc

    final_size = part.stat().st_size
    if final_size != expected_size:
        raise SystemExit(
            f"download incomplete: got {final_size}, expected {expected_size}. "
            f"Keep {part} and rerun to resume."
        )

    part.replace(dest)
    print(f"done: {dest} ({fmt_size(final_size)})")


def main() -> int:
    parser = argparse.ArgumentParser(description="Download SRR32296234 SRA fixture.")
    parser.add_argument("--dest", default=str(DEFAULT_DEST), help="destination path")
    parser.add_argument("--url", default=URL, help="source URL")
    parser.add_argument("--force", action="store_true", help="replace existing destination")
    parser.add_argument(
        "--skip-head",
        action="store_true",
        help="skip HEAD size check before downloading",
    )
    args = parser.parse_args()

    expected_size = EXPECTED_SIZE
    if not args.skip_head:
        size = remote_size(args.url)
        if size != EXPECTED_SIZE:
            raise SystemExit(f"remote size changed: got {size}, expected {EXPECTED_SIZE}")
        expected_size = size

    download(args.url, Path(args.dest).resolve(), expected_size, args.force)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
