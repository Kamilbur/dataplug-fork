"""Tests for the SRA → FASTQ dataplug backend.

S3 is mocked via moto — no MinIO required for these tests.
VDB reads SRA files by accession name from the current working directory,
so each test runs with CWD set to SRA_FIXTURES_DIR.

Output format of SRASlice.lazy_get():
    Yields one tuple[str, ...] per spot. Each element of the tuple is a
    complete FASTQ quadruple string (4 lines joined by '\\n', with a trailing
    '\\n'), corresponding to one split read. For non-split mode (split=False)
    the tuple always has exactly one element.

Reference hashes in fixtures/hashes/ were produced by fasterq-dump with
HASH_STEP=1000 records per SHA-256 segment:
    splits[0]  — hashes for _1.fastq  (read 1 from --split-3)
    splits[1]  — hashes for _2.fastq  (read 2 from --split-3)
    concats[0] — hashes for .fastq    (--concatenate-reads)
"""

import hashlib

import boto3
import pytest
from moto import mock_aws

from conftest import (
    HASH_STEP,
    META_BUCKET,
    S3_BUCKET,
    SRA_FIXTURES_DIR,
    load_reference_hashes,
    sra_accessions_with_fixtures,
)

# ---------------------------------------------------------------------------
# Skip entire module if native vdb library is not available
# ---------------------------------------------------------------------------

try:
    from dataplug import CloudObject
    from dataplug.formats.genomics.sra import SRA, partition_chunks_strategy
    VDB_AVAILABLE = True
except (ImportError, RuntimeError) as e:
    import traceback
    print(e)
    traceback.print_exc()
    VDB_AVAILABLE = False

pytestmark = pytest.mark.skipif(not VDB_AVAILABLE, reason="vdb native library not available")

CHUNK_COUNTS = [1, 5, 10, 20, 27]


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(autouse=True)
def run_from_sra_fixtures_dir(monkeypatch):
    """VDB opens SRA files by accession name relative to CWD."""
    monkeypatch.chdir(SRA_FIXTURES_DIR)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _upload_sra_to_mock_s3(acc: str) -> str:
    """Create buckets, upload local SRA fixture, return S3 URI."""
    s3 = boto3.client("s3", region_name="us-east-1")
    for bucket in (S3_BUCKET, META_BUCKET):
        s3.create_bucket(Bucket=bucket)
    s3.upload_file(str(SRA_FIXTURES_DIR / acc), S3_BUCKET, acc)
    return f"s3://{S3_BUCKET}/{acc}"


def _preprocess_and_partition(uri: str, num_chunks: int):
    co = CloudObject.from_s3(SRA, uri, metadata_bucket=META_BUCKET)
    co.preprocess()
    return co.partition(partition_chunks_strategy, num_chunks=num_chunks)


def _collect_split_records(chunks) -> tuple[list[str], list[str]]:
    """Collect per-read FASTQ records from split mode.

    Returns (read1_records, read2_records). Each record is a complete
    FASTQ quadruple string with trailing newline, as yielded by
    FastqLineQuadruple.__str__().
    """
    read1, read2 = [], []
    for chunk in chunks:
        for reads in chunk.lazy_get(split=True):
            if len(reads) > 0:
                read1.append(reads[0])
            if len(reads) > 1:
                read2.append(reads[1])
    return read1, read2


def _collect_concat_records(chunks) -> list[str]:
    """Collect FASTQ records from concat mode (split=False)."""
    records = []
    for chunk in chunks:
        for reads in chunk.lazy_get(split=False):
            records.extend(reads)
    return records


def _hash_records(records: list[str], step: int = HASH_STEP) -> list[str]:
    """SHA-256 over batches of `step` records, matching prepare_tests generation."""
    hashes = []
    for start in range(0, len(records), step):
        sha = hashlib.sha256()
        for record in records[start : start + step]:
            sha.update(record.encode("utf-8"))
        hashes.append(sha.hexdigest())
    return hashes


@pytest.mark.parametrize("acc", sra_accessions_with_fixtures())
@mock_aws
def test_output_is_chunk_count_invariant(acc):
    """Concatenated output across all chunks must be identical regardless of chunk count."""
    uri = _upload_sra_to_mock_s3(acc)
    co = CloudObject.from_s3(SRA, uri, metadata_bucket=META_BUCKET)
    co.preprocess()

    all_records = {
        nc: _collect_concat_records(co.partition(partition_chunks_strategy, num_chunks=nc))
        for nc in CHUNK_COUNTS
    }

    baseline_nc = CHUNK_COUNTS[0]
    baseline = all_records[baseline_nc]
    for nc, records in all_records.items():
        assert records == baseline, (
            f"num_chunks={nc} produced {len(records)} records, "
            f"expected {len(baseline)} (same as num_chunks={baseline_nc})"
        )


@pytest.mark.parametrize("acc", sra_accessions_with_fixtures())
@mock_aws
def test_fastq_structure_is_valid(acc):
    """Every output record must be a valid FASTQ quadruple."""
    uri = _upload_sra_to_mock_s3(acc)
    chunks = _preprocess_and_partition(uri, num_chunks=5)

    for chunk in chunks:
        for reads in chunk.lazy_get(split=True):
            for record in reads:
                lines = record.rstrip("\n").split("\n")
                assert len(lines) == 4, f"Expected 4 lines per record, got {len(lines)}"
                header1, seq, header2, qual = lines
                assert header1.startswith("@"), f"Line 0 must start with '@': {header1!r}"
                assert header2.startswith("+"), f"Line 2 must start with '+': {header2!r}"
                assert len(seq) == len(qual), (
                    f"Seq length ({len(seq)}) != qual length ({len(qual)})"
                )

