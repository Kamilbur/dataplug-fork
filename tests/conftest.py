import json
import os
import shutil
import urllib.request
from pathlib import Path

import boto3
import pytest
from moto import mock_aws

FIXTURES_DIR = Path(__file__).parent / "fixtures"
SRA_FIXTURES_DIR = FIXTURES_DIR / "sra"
HASHES_DIR = FIXTURES_DIR / "hashes"
FIXTURE_MANIFEST = Path(__file__).parent / "prepare_tests" / "test-files.json"

S3_ENDPOINT = os.getenv("S3_ENDPOINT", "http://localhost:9000")
S3_ACCESS_KEY = os.getenv("S3_ACCESS_KEY", "minioadmin")
S3_SECRET_KEY = os.getenv("S3_SECRET_KEY", "minioadmin")
S3_BUCKET = os.getenv("S3_BUCKET", "test-sra")
META_BUCKET = S3_BUCKET + ".meta"


_DOWNLOAD_CHUNK = 256 * 1024 * 1024  # 256 MB


def pytest_addoption(parser):
    parser.addoption(
        "--download-fixtures",
        action="store_true",
        default=False,
        help=(
            "Download missing SRA fixture files from NCBI before running tests. "
            "Files are saved to tests/fixtures/fixtures/."
        ),
    )


def load_reference_hashes(acc: str) -> dict:
    path = HASHES_DIR / f"{acc}.json"
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def sra_accessions_with_fixtures() -> list[str]:
    """Return accessions that have a local SRA binary in fixtures/fixtures/."""
    if not SRA_FIXTURES_DIR.exists():
        return []
    return sorted(p.name for p in SRA_FIXTURES_DIR.iterdir() if p.is_file())


def _load_manifest() -> list[dict]:
    return json.loads(FIXTURE_MANIFEST.read_text(encoding="utf-8"))["files"]


def _download_sra_file(acc: str, url: str) -> None:
    SRA_FIXTURES_DIR.mkdir(parents=True, exist_ok=True)
    dest = SRA_FIXTURES_DIR / acc
    print(f"\nDownloading {acc} from {url} ...")
    with urllib.request.urlopen(url) as resp, dest.open("wb") as f:
        shutil.copyfileobj(resp, f, length=_DOWNLOAD_CHUNK)
    print(f"Saved {acc} ({dest.stat().st_size // 1024 // 1024} MB)")


@pytest.fixture
def s3_uri(request):
    """Upload local SRA fixture to mocked S3, yield S3 URI. request.param = accession."""
    acc = request.param
    with mock_aws():
        s3 = boto3.client("s3", region_name="us-east-1")
        for bucket in (S3_BUCKET, META_BUCKET):
            s3.create_bucket(Bucket=bucket)
        s3.upload_file(str(SRA_FIXTURES_DIR / acc), S3_BUCKET, acc)
        yield f"s3://{S3_BUCKET}/{acc}"


@pytest.fixture(scope="session", autouse=True)
def sra_fixtures(pytestconfig):
    """Ensure local SRA fixture files are present.

    Without --download-fixtures: uses whatever is already in fixtures/fixtures/.
    With --download-fixtures: downloads any missing files listed in the manifest.
    """
    if not pytestconfig.getoption("--download-fixtures"):
        return

    manifest = {entry["acc"]: entry["sra-url"] for entry in _load_manifest()}
    expected = set(manifest) & {p.stem for p in HASHES_DIR.glob("*.json")}
    missing = {acc for acc in expected if not (SRA_FIXTURES_DIR / acc).exists()}

    for acc in sorted(missing):
        _download_sra_file(acc, manifest[acc])
