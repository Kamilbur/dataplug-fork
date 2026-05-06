"""Tests for the MGB -> FASTQ dataplug backend.

S3 is mocked via moto. GENIE reads .mgb files by local filename through
the dataplug shim, so native libgenie is loaded only inside active tests.
"""

from __future__ import annotations

import os
from pathlib import Path

import boto3
import pytest
from moto import mock_aws

from conftest import FIXTURES_DIR

try:
    from dataplug import CloudObject
    from dataplug.formats.genomics.mgb import MGB, partition_chunks_strategy

    MGB_FORMAT_AVAILABLE = True
except (ImportError, RuntimeError):
    MGB_FORMAT_AVAILABLE = False


MGB_BUCKET = os.getenv("MGB_S3_BUCKET", "test-mgb")
MGB_META_BUCKET = MGB_BUCKET + ".meta"
MGB_FIXTURES_DIR = Path(os.getenv("MGB_FIXTURES_DIR", FIXTURES_DIR / "mgb"))
FASTQ_FIXTURES_DIR = Path(os.getenv("FASTQ_FIXTURES_DIR", FIXTURES_DIR / "fastq"))
GENIE_CLONE_FIXTURES_DIR = (
    Path(__file__).resolve().parents[1]
    / ".clones"
    / "genie"
    / "data"
    / "mpeg-g"
    / "conformance"
)
CHUNK_COUNTS = [1, 2, 3, 5]


def _genie_library_exists() -> bool:
    env_path = os.environ.get("GENIE_SO_PATH")
    if env_path:
        return Path(env_path).exists()

    repo_root = Path(__file__).resolve().parents[1]
    candidates = [
        repo_root / "dataplug" / "formats" / "genomics" / "mgb" / "internals" / "libgenie.so",
        repo_root / ".clones" / "genie" / "build" / "lib" / "libgenie.so",
    ]
    return any(path.exists() for path in candidates)


def mgb_files_with_fixtures() -> list[Path]:
    paths = []
    if MGB_FIXTURES_DIR.exists():
        paths.extend(sorted(MGB_FIXTURES_DIR.glob("*.mgb")))
    if not paths and GENIE_CLONE_FIXTURES_DIR.exists():
        paths.extend(sorted(GENIE_CLONE_FIXTURES_DIR.glob("*/*.mgb")))
    return paths


def mgb_fastq_fixture_pairs() -> list[tuple[Path, Path]]:
    """Return MGB fixtures that have an accession-matched FASTQ fixture."""
    pairs = []
    for mgb_path in mgb_files_with_fixtures():
        fastq_path = FASTQ_FIXTURES_DIR / f"{mgb_path.stem}.fastq"
        if fastq_path.exists():
            pairs.append((mgb_path, fastq_path))
    return pairs


MGB_FASTQ_FIXTURE_PAIRS = mgb_fastq_fixture_pairs()


MGB_AVAILABLE = MGB_FORMAT_AVAILABLE and _genie_library_exists()

pytestmark = pytest.mark.skipif(
    not MGB_AVAILABLE,
    reason="GENIE native library not available",
)


@pytest.fixture
def mgb_s3_uri(request):
    """Upload local MGB fixture to mocked S3, yield S3 URI."""
    fixture = Path(request.param)
    with mock_aws():
        s3 = boto3.client("s3", region_name="us-east-1")
        for bucket in (MGB_BUCKET, MGB_META_BUCKET):
            s3.create_bucket(Bucket=bucket)
        s3.upload_file(str(fixture), MGB_BUCKET, fixture.name)
        yield f"s3://{MGB_BUCKET}/{fixture.name}"


def _preprocess_and_partition(uri: str, num_chunks: int):
    co = CloudObject.from_s3(MGB, uri, metadata_bucket=MGB_META_BUCKET)
    co.preprocess()
    return co.partition(partition_chunks_strategy, num_chunks=num_chunks)


def _collect_records(chunks) -> list[str]:
    records = []
    for chunk in chunks:
        for reads in chunk.lazy_get():
            records.extend(reads)
    return records


def _load_fastq_records(path: Path) -> list[str]:
    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    assert len(lines) % 4 == 0, f"Expected FASTQ line count to be divisible by 4: {path}"
    return ["".join(lines[idx:idx + 4]) for idx in range(0, len(lines), 4)]


def _assert_records_match(actual: list[str], expected: list[str]) -> None:
    assert len(actual) == len(expected), (
        f"Record count mismatch: got {len(actual)}, expected {len(expected)}"
    )
    for idx, (actual_record, expected_record) in enumerate(zip(actual, expected)):
        if actual_record != expected_record:
            actual_header = actual_record.split("\n", 1)[0]
            expected_header = expected_record.split("\n", 1)[0]
            assert actual_record == expected_record, (
                f"Record {idx} mismatch: got {actual_header!r}, "
                f"expected {expected_header!r}"
            )


@pytest.mark.parametrize("mgb_s3_uri", mgb_files_with_fixtures(), indirect=True)
def test_output_is_chunk_count_invariant(mgb_s3_uri):
    """Output across all chunks must be identical regardless of chunk count."""
    co = CloudObject.from_s3(MGB, mgb_s3_uri, metadata_bucket=MGB_META_BUCKET)
    co.preprocess()

    all_records = {
        chunk_count: _collect_records(
            co.partition(partition_chunks_strategy, num_chunks=chunk_count)
        )
        for chunk_count in CHUNK_COUNTS
    }

    baseline_count = CHUNK_COUNTS[0]
    baseline = all_records[baseline_count]
    for chunk_count, records in all_records.items():
        assert records == baseline, (
            f"num_chunks={chunk_count} produced {len(records)} records, "
            f"expected {len(baseline)} (same as num_chunks={baseline_count})"
        )


@pytest.mark.parametrize("mgb_s3_uri", mgb_files_with_fixtures(), indirect=True)
def test_fastq_structure_is_valid(mgb_s3_uri):
    """Every output record must be a valid FASTQ quadruple."""
    chunks = _preprocess_and_partition(mgb_s3_uri, num_chunks=5)

    saw_records = False
    for chunk in chunks:
        for reads in chunk.lazy_get():
            for record in reads:
                saw_records = True
                lines = record.rstrip("\n").split("\n")
                assert len(lines) == 4, f"Expected 4 lines per record, got {len(lines)}"
                header1, seq, header2, qual = lines
                assert header1.startswith("@"), f"Line 0 must start with '@': {header1!r}"
                assert header2.startswith("+"), f"Line 2 must start with '+': {header2!r}"
                assert len(seq) == len(qual), (
                    f"Seq length ({len(seq)}) != qual length ({len(qual)})"
                )

    assert saw_records, "Expected at least one FASTQ record"


@pytest.mark.parametrize(
    ("mgb_s3_uri", "expected_fastq_path"),
    MGB_FASTQ_FIXTURE_PAIRS,
    indirect=["mgb_s3_uri"],
    ids=[path.stem for path, _ in MGB_FASTQ_FIXTURE_PAIRS],
)
def test_output_matches_fastq_fixture(mgb_s3_uri, expected_fastq_path):
    """MGB output must match accession-matched FASTQ fixtures exactly."""
    chunks = _preprocess_and_partition(mgb_s3_uri, num_chunks=5)

    _assert_records_match(
        _collect_records(chunks),
        _load_fastq_records(expected_fastq_path),
    )
