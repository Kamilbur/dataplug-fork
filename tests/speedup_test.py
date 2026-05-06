"""Measure SRA preprocessing speedup across download worker counts.

This is a local benchmark helper, not a pytest test.

Examples:
    source venv/bin/activate
    python3 tests/speedup_test.py --s3-uri s3://test-sra/SRR19392985
    python3 tests/speedup_test.py --s3-uri s3://test-sra/SRR19392985 --workers 1,2,3,5 --repeats 2
    python3 tests/speedup_test.py --s3-uri s3://test-sra/SRR19392985 --s3-config-json '{"endpoint_url":"http://localhost:9000","region_name":"us-east-1","credentials":{"AccessKeyId":"minioadmin","SecretAccessKey":"minioadmin"}}'

Outputs:
    benchmark_results/sra_speedup_raw.csv
    benchmark_results/sra_speedup_raw.jsonl
    benchmark_results/sra_speedup_summary.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import selectors
import subprocess
import sys
import time
from datetime import UTC, datetime
from pathlib import Path
from statistics import mean, stdev

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = REPO_ROOT / "benchmark_results"
RAW_FIELDS = [
    "timestamp_utc",
    "accession",
    "s3_uri",
    "size_bytes",
    "cpu_count",
    "workers",
    "walk_workers",
    "repeat",
    "step",
    "elapsed_s",
    "maxrss_mb",
    "total_lines",
    "pread_groups",
    "pread_intervals",
    "mmap_count",
    "status",
    "error",
]
SUMMARY_FIELDS = [
    "accession",
    "s3_uri",
    "size_bytes",
    "cpu_count",
    "workers",
    "walk_workers",
    "runs",
    "step",
    "mean_elapsed_s",
    "stdev_elapsed_s",
    "min_elapsed_s",
    "max_elapsed_s",
    "mean_maxrss_mb",
    "speedup_vs_1_worker",
]


def parse_workers(value: str | None, cpu_count: int) -> list[int]:
    if value:
        workers = sorted({int(part) for part in value.split(",") if part.strip()})
    else:
        # Default never uses all visible CPUs. User can request it explicitly with
        # --include-max or --workers.
        workers = list(range(1, max(1, cpu_count)))
    if any(worker < 1 for worker in workers):
        raise argparse.ArgumentTypeError("workers must be positive integers")
    return workers


def append_csv(path: Path, fieldnames: list[str], row: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    exists = path.exists()
    with path.open("a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if not exists:
            writer.writeheader()
        writer.writerow({key: row.get(key, "") for key in fieldnames})


def write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def append_jsonl(path: Path, row: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as f:
        f.write(json.dumps(row, sort_keys=True) + "\n")


def load_s3_config(args: argparse.Namespace) -> dict:
    if args.s3_config_json:
        return json.loads(args.s3_config_json)
    if args.s3_config_file:
        return json.loads(Path(args.s3_config_file).read_text(encoding="utf-8"))

    config: dict = {}
    if args.s3_endpoint:
        config["endpoint_url"] = args.s3_endpoint
    if args.s3_region:
        config["region_name"] = args.s3_region
    if bool(args.s3_access_key) != bool(args.s3_secret_key):
        raise ValueError("--s3-access-key and --s3-secret-key must be provided together")
    if args.s3_access_key and args.s3_secret_key:
        credentials = {
            "AccessKeyId": args.s3_access_key,
            "SecretAccessKey": args.s3_secret_key,
        }
        if args.s3_session_token:
            credentials["SessionToken"] = args.s3_session_token
        config["credentials"] = credentials
    if args.s3_unsigned:
        config.setdefault("botocore_config_kwargs", {})["signature_version"] = "unsigned"
    return config


def child_run(args: argparse.Namespace) -> int:
    sys.path.insert(0, str(REPO_ROOT))

    from dataplug import CloudObject
    import dataplug.formats.genomics.sra.experimental as exp
    from dataplug.formats.genomics.sra.internals.download import (
        parallel_download_to_buffer as real_download,
    )

    s3_config = load_s3_config(args)
    accession = args.accession or args.s3_uri.rstrip("/").split("/")[-1].split(".")[0]
    worker_count = int(args.workers)

    def download_with_worker_count(cloud_object):
        return real_download(cloud_object, max_workers=worker_count)

    exp.parallel_download_to_buffer = download_with_worker_count

    co = CloudObject.from_s3(
        exp.SRA,
        args.s3_uri,
        metadata_bucket=args.metadata_bucket,
        s3_config=s3_config,
    )
    print(
        f"child start workers={worker_count} walk_workers={args.walk_workers} "
        f"repeat={args.repeat} size_bytes={co.size} step={args.step}",
        flush=True,
    )
    started = time.perf_counter()
    md = exp.preprocess_sra(
        co,
        step=args.step,
        walk_workers=args.walk_workers,
        walk_start_method=args.walk_start_method,
        walk_overlap_steps=args.walk_overlap_steps,
    )
    elapsed = time.perf_counter() - started
    print(f"child done workers={worker_count} repeat={args.repeat} elapsed_s={elapsed:.6f}", flush=True)

    try:
        import resource

        maxrss_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    except Exception:
        maxrss_mb = 0.0

    attrs = md.attributes
    total_lines = int(attrs["total_lines"])
    if args.expected_lines is not None and total_lines != args.expected_lines:
        raise RuntimeError(
            f"total_lines mismatch: expected {args.expected_lines}, got {total_lines}"
        )

    row = {
        "timestamp_utc": datetime.now(UTC).isoformat(timespec="seconds"),
        "accession": accession,
        "s3_uri": args.s3_uri,
        "size_bytes": co.size,
        "cpu_count": os.cpu_count() or 1,
        "workers": worker_count,
        "walk_workers": args.walk_workers,
        "repeat": args.repeat,
        "step": args.step,
        "elapsed_s": round(elapsed, 6),
        "maxrss_mb": round(maxrss_mb, 3),
        "total_lines": total_lines,
        "pread_groups": len(attrs["preads"]),
        "pread_intervals": sum(len(value) for value in attrs["preads"].values()),
        "mmap_count": len(attrs["mmaps"]),
        "status": "ok",
        "error": "",
    }
    print(json.dumps(row, sort_keys=True), flush=True)
    return 0


def run_child(args: argparse.Namespace, workers: int, repeat: int) -> dict:
    cmd = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--child-run",
        "--s3-uri",
        args.s3_uri,
        "--accession",
        args.accession,
        "--metadata-bucket",
        args.metadata_bucket,
        "--workers",
        str(workers),
        "--repeat",
        str(repeat),
        "--step",
        str(args.step),
    ]
    cmd.extend(["--walk-workers", str(args.walk_workers)])
    cmd.extend(["--walk-start-method", args.walk_start_method])
    cmd.extend(["--walk-overlap-steps", str(args.walk_overlap_steps)])
    if args.s3_config_json:
        cmd.extend(["--s3-config-json", args.s3_config_json])
    if args.s3_config_file:
        cmd.extend(["--s3-config-file", args.s3_config_file])
    if args.s3_endpoint:
        cmd.extend(["--s3-endpoint", args.s3_endpoint])
    if args.s3_region:
        cmd.extend(["--s3-region", args.s3_region])
    if args.s3_access_key:
        cmd.extend(["--s3-access-key", args.s3_access_key])
    if args.s3_secret_key:
        cmd.extend(["--s3-secret-key", args.s3_secret_key])
    if args.s3_session_token:
        cmd.extend(["--s3-session-token", args.s3_session_token])
    if args.s3_unsigned:
        cmd.append("--s3-unsigned")
    if args.expected_lines is not None:
        cmd.extend(["--expected-lines", str(args.expected_lines)])

    proc = subprocess.Popen(
        cmd,
        cwd=REPO_ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        bufsize=1,
    )
    assert proc.stdout is not None
    assert proc.stderr is not None

    stdout_lines: list[str] = []
    stderr_lines: list[str] = []
    json_row: dict | None = None
    selector = selectors.DefaultSelector()
    selector.register(proc.stdout, selectors.EVENT_READ, "stdout")
    selector.register(proc.stderr, selectors.EVENT_READ, "stderr")
    last_heartbeat = time.monotonic()

    while selector.get_map():
        events = selector.select(timeout=1.0)
        if not events:
            if proc.poll() is None and time.monotonic() - last_heartbeat >= args.heartbeat_interval:
                print(f"still running workers={workers} repeat={repeat}", flush=True)
                last_heartbeat = time.monotonic()
            continue
        for key, _ in events:
            line = key.fileobj.readline()
            if line == "":
                selector.unregister(key.fileobj)
                continue
            stream_name = key.data
            clean = line.rstrip("\n")
            if stream_name == "stdout":
                stdout_lines.append(clean)
                if clean.startswith("{") and clean.endswith("}"):
                    json_row = json.loads(clean)
                else:
                    print(f"[child stdout] {clean}", flush=True)
            else:
                stderr_lines.append(clean)
                print(f"[child stderr] {clean}", flush=True)

    returncode = proc.wait()
    if returncode == 0 and json_row is not None:
        return json_row

    return {
        "timestamp_utc": datetime.now(UTC).isoformat(timespec="seconds"),
        "accession": args.accession,
        "s3_uri": args.s3_uri,
        "size_bytes": "",
        "cpu_count": os.cpu_count() or 1,
        "workers": workers,
        "walk_workers": args.walk_workers,
        "repeat": repeat,
        "step": args.step,
        "elapsed_s": "",
        "maxrss_mb": "",
        "total_lines": "",
        "pread_groups": "",
        "pread_intervals": "",
        "mmap_count": "",
        "status": f"failed:{returncode}",
        "error": ("\n".join(stderr_lines).strip() or "\n".join(stdout_lines).strip())[-2000:],
    }


def summarize(rows: list[dict]) -> list[dict]:
    ok_rows = [row for row in rows if row.get("status") == "ok"]
    by_workers: dict[int, list[dict]] = {}
    for row in ok_rows:
        by_workers.setdefault(int(row["workers"]), []).append(row)

    baseline_rows = by_workers.get(1, [])
    baseline = mean(float(row["elapsed_s"]) for row in baseline_rows) if baseline_rows else None

    summary = []
    for workers in sorted(by_workers):
        group = by_workers[workers]
        elapsed = [float(row["elapsed_s"]) for row in group]
        maxrss = [float(row["maxrss_mb"]) for row in group]
        row0 = group[0]
        summary.append(
            {
                "accession": row0["accession"],
                "s3_uri": row0["s3_uri"],
                "size_bytes": row0["size_bytes"],
                "cpu_count": row0["cpu_count"],
                "workers": workers,
                "walk_workers": row0["walk_workers"],
                "runs": len(group),
                "step": row0["step"],
                "mean_elapsed_s": round(mean(elapsed), 6),
                "stdev_elapsed_s": round(stdev(elapsed), 6) if len(elapsed) > 1 else 0.0,
                "min_elapsed_s": round(min(elapsed), 6),
                "max_elapsed_s": round(max(elapsed), 6),
                "mean_maxrss_mb": round(mean(maxrss), 3),
                "speedup_vs_1_worker": (
                    round(baseline / mean(elapsed), 6) if baseline else ""
                ),
            }
        )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Benchmark SRA preprocess_sra speedup across worker counts."
    )
    parser.add_argument("--s3-uri", required=True, help="SRA object URI, e.g. s3://bucket/SRR19392985")
    parser.add_argument("--accession", default="", help="Accession name; default: S3 key basename")
    parser.add_argument("--metadata-bucket", default="", help="Metadata bucket; default: data bucket + .meta")
    parser.add_argument("--s3-config-json", default="", help="JSON dict passed as CloudObject s3_config")
    parser.add_argument("--s3-config-file", default="", help="Path to JSON dict passed as CloudObject s3_config")
    parser.add_argument("--s3-endpoint", default=os.getenv("S3_ENDPOINT", ""))
    parser.add_argument("--s3-region", default=os.getenv("AWS_DEFAULT_REGION", os.getenv("AWS_REGION", "")))
    parser.add_argument("--s3-access-key", default=os.getenv("S3_ACCESS_KEY", os.getenv("AWS_ACCESS_KEY_ID", "")))
    parser.add_argument("--s3-secret-key", default=os.getenv("S3_SECRET_KEY", os.getenv("AWS_SECRET_ACCESS_KEY", "")))
    parser.add_argument("--s3-session-token", default=os.getenv("AWS_SESSION_TOKEN", ""))
    parser.add_argument("--s3-unsigned", action="store_true", help="Use unsigned S3 requests")
    parser.add_argument("--workers", default="", help="Comma list, e.g. 1,2,3,5")
    parser.add_argument("--walk-workers", type=int, default=1, help="Processes for VDB row-walk phase")
    parser.add_argument("--walk-start-method", default="spawn", choices=("spawn", "fork", "forkserver"))
    parser.add_argument("--walk-overlap-steps", type=int, default=0, help="Overlapped sampled rows per worker boundary")
    parser.add_argument("--include-max", action="store_true", help="Include os.cpu_count() in default workers")
    parser.add_argument("--repeats", type=int, default=1, help="Runs per worker count")
    parser.add_argument("--step", type=int, default=250, help="preprocess_sra sampling step")
    parser.add_argument("--expected-lines", type=int, default=None)
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    parser.add_argument("--heartbeat-interval", type=float, default=30.0, help="Seconds between parent progress pings")
    parser.add_argument("--child-run", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--repeat", type=int, default=1, help=argparse.SUPPRESS)
    args = parser.parse_args()

    if not args.accession:
        args.accession = args.s3_uri.rstrip("/").split("/")[-1].split(".")[0]
    if not args.metadata_bucket:
        bucket = args.s3_uri.removeprefix("s3://").split("/", 1)[0]
        args.metadata_bucket = bucket + ".meta"
    if args.repeats < 1:
        parser.error("--repeats must be >= 1")
    if args.step < 1:
        parser.error("--step must be >= 1")
    if args.walk_workers < 1:
        parser.error("--walk-workers must be >= 1")
    if args.walk_overlap_steps < 0:
        parser.error("--walk-overlap-steps must be >= 0")

    if args.child_run:
        return child_run(args)

    try:
        load_s3_config(args)
    except ValueError as exc:
        parser.error(str(exc))

    cpu_count = os.cpu_count() or 1
    workers = parse_workers(args.workers, cpu_count)
    if not args.workers and args.include_max and cpu_count not in workers:
        workers.append(cpu_count)
    workers = sorted(workers)

    output_dir = Path(args.output_dir).resolve()
    raw_csv = output_dir / "sra_speedup_raw.csv"
    raw_jsonl = output_dir / "sra_speedup_raw.jsonl"
    summary_csv = output_dir / "sra_speedup_summary.csv"

    print(f"s3_uri={args.s3_uri}")
    print(f"accession={args.accession}")
    print(f"metadata_bucket={args.metadata_bucket}")
    print(f"cpu_count={cpu_count}")
    print(f"workers={workers}")
    print(f"walk_workers={args.walk_workers}")
    print(f"walk_start_method={args.walk_start_method}")
    print(f"walk_overlap_steps={args.walk_overlap_steps}")
    print(f"repeats={args.repeats}")
    print(f"raw_csv={raw_csv}")
    print(f"raw_jsonl={raw_jsonl}")
    print(f"summary_csv={summary_csv}")

    rows = []
    for repeat in range(1, args.repeats + 1):
        for worker_count in workers:
            print(f"run workers={worker_count} repeat={repeat}", flush=True)
            row = run_child(args, worker_count, repeat)
            rows.append(row)
            append_csv(raw_csv, RAW_FIELDS, row)
            append_jsonl(raw_jsonl, row)
            if row["status"] == "ok":
                print(
                    "ok "
                    f"workers={worker_count} "
                    f"elapsed_s={row['elapsed_s']} "
                    f"maxrss_mb={row['maxrss_mb']}"
                )
            else:
                print(f"fail workers={worker_count} status={row['status']}")

    summary = summarize(rows)
    write_csv(summary_csv, SUMMARY_FIELDS, summary)

    print("\nsummary")
    print("workers,mean_elapsed_s,speedup_vs_1_worker,mean_maxrss_mb")
    for row in summary:
        print(
            f"{row['workers']},{row['mean_elapsed_s']},"
            f"{row['speedup_vs_1_worker']},{row['mean_maxrss_mb']}"
        )
    return 0 if all(row["status"] == "ok" for row in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
