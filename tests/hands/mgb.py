"""Manual run of MGB FASTQ validation for SRR19392985 for fast prototyping.

Usage:
    python3 tests/hands/mgb.py
"""

import sys
from pathlib import Path

import boto3
from moto import mock_aws


def main():
    bucket = "test-mgb"
    meta_bucket = bucket + ".meta"
    accession = "SRR19392985"
    filename = f"{accession}.mgb"
    uri = f"s3://{bucket}/{filename}"
    fixture = Path(__file__).parents[1] / "fixtures" / "mgb" / filename

    with mock_aws():
        s3 = boto3.client("s3", region_name="us-east-1")
        for bucket_name in (bucket, meta_bucket):
            s3.create_bucket(Bucket=bucket_name)
        s3.upload_file(str(fixture), bucket, filename)

        from dataplug import CloudObject
        from dataplug.formats.genomics.mgb import (
            MGB,
            partition_chunks_strategy,
        )

        co = CloudObject.from_s3(MGB, uri, metadata_bucket=meta_bucket)
        print("Preprocessing...")
        co.preprocess(force=True)
        chunks = co.partition(partition_chunks_strategy, num_chunks=5)
        print(f"Partitioned into {len(chunks)} chunks")

        total_records = 0
        errors = []

        empty_reads = False
        for chunk_idx, chunk in enumerate(chunks):
            for reads in chunk.lazy_get():
                if len(reads) == 0:
                    empty_reads = True
                for record in reads:
                    lines = record.rstrip("\n").split("\n")
                    print(lines)
                    total_records += 1
                    if len(lines) != 4:
                        errors.append(f"chunk {chunk_idx}: expected 4 lines, got {len(lines)}")
                        continue
                    header1, seq, header2, qual = lines
                    if not header1.startswith("@"):
                        errors.append(f"chunk {chunk_idx}: line 0 must start with '@': {header1!r}")
                    if not header2.startswith("+"):
                        errors.append(f"chunk {chunk_idx}: line 2 must start with '+': {header2!r}")
                    if len(seq) != len(qual):
                        errors.append(
                            f"chunk {chunk_idx}: seq len {len(seq)} != qual len {len(qual)}"
                        )

        if empty_reads:
            errors.append("empty_reads")
        print(f"Validated {total_records} records across 5 chunks")

        if errors:
            print(f"FAIL — {len(errors)} error(s):")
            for e in errors:
                print(f"  {e}")
            sys.exit(1)
        else:
            print("PASS")


if __name__ == "__main__":
    main()
