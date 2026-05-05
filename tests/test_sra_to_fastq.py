"""Tests for the SRA → FASTQ dataplug backend.

S3 is mocked via moto. VDB reads SRA files by accession name from
the current working directory, so each test runs with CWD set
to SRA_FIXTURES_DIR.

"""

import pytest

from conftest import (
    META_BUCKET,
    sra_accessions_with_fixtures,
)

try:
    from dataplug import CloudObject
    from dataplug.formats.genomics.sra.experimental import (
        SRA,
        partition_chunks_strategy,
    )

    SRA_AVAILABLE = True
except (ImportError, RuntimeError):
    SRA_AVAILABLE = False

pytestmark = pytest.mark.skipif(not SRA_AVAILABLE, reason="vdb native library not available")

CHUNK_COUNTS = [1, 5, 10, 20, 27]


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


@pytest.mark.parametrize("s3_uri", sra_accessions_with_fixtures(), indirect=True)
def test_output_is_chunk_count_invariant(s3_uri):
    """Concatenated output across all chunks must be identical regardless of chunk count."""
    co = CloudObject.from_s3(SRA, s3_uri, metadata_bucket=META_BUCKET)
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


@pytest.mark.parametrize("s3_uri", sra_accessions_with_fixtures(), indirect=True)
def test_fastq_structure_is_valid(s3_uri):
    """Every output record must be a valid FASTQ quadruple."""
    chunks = _preprocess_and_partition(s3_uri, num_chunks=5)

    for chunk in chunks:
        for reads in chunk.lazy_get():
            for record in reads:
                lines = record.rstrip("\n").split("\n")
                assert len(lines) == 4, f"Expected 4 lines per record, got {len(lines)}"
                header1, seq, header2, qual = lines
                assert header1.startswith("@"), f"Line 0 must start with '@': {header1!r}"
                assert header2.startswith("+"), f"Line 2 must start with '+': {header2!r}"
                assert len(seq) == len(qual), (
                    f"Seq length ({len(seq)}) != qual length ({len(qual)})"
                )
