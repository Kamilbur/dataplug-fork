"""Local-file-backed S3-compatible storage stub for tests.

Avoids moto's in-memory blowup when the SRA fixture is multi-GB.
Implements just enough of the boto3 S3 client interface for
CloudObject + experimental preprocess + partition.
"""

from __future__ import annotations

import io
import os
import pickle
from threading import Lock

import botocore.exceptions


def _client_error(code: str, op: str) -> botocore.exceptions.ClientError:
    return botocore.exceptions.ClientError(
        {"Error": {"Code": code, "Message": code}, "ResponseMetadata": {"HTTPStatusCode": 404}},
        op,
    )


class LocalFileS3:
    """A minimal S3 client backed by a local file for the data bucket and an
    in-process dict for the meta bucket(s).
    """

    def __init__(self, data_path: str, data_bucket: str, data_key: str):
        self._data_path = str(data_path)
        self._data_bucket = data_bucket
        self._data_key = data_key
        self._data_size = os.path.getsize(self._data_path)
        self._meta: dict[str, dict[str, bytes]] = {}
        self._buckets: set[str] = {data_bucket}
        self._lock = Lock()

    def _read_data_range(self, left: int, right: int | None) -> bytes:
        with open(self._data_path, "rb") as f:
            f.seek(left)
            if right is None:
                return f.read()
            return f.read(right - left + 1)

    def _is_data(self, bucket: str, key: str) -> bool:
        return bucket == self._data_bucket and key == self._data_key

    def head_object(self, Bucket, Key, **_):
        if self._is_data(Bucket, Key):
            return {
                "ContentLength": self._data_size,
                "Metadata": {},
                "ResponseMetadata": {"HTTPStatusCode": 200},
            }
        with self._lock:
            data = self._meta.get(Bucket, {}).get(Key)
        if data is None:
            raise _client_error("404", "HeadObject")
        return {
            "ContentLength": len(data),
            "Metadata": {},
            "ResponseMetadata": {"HTTPStatusCode": 200},
        }

    def head_bucket(self, Bucket, **_):
        with self._lock:
            if Bucket not in self._buckets:
                raise _client_error("404", "HeadBucket")
        return {"ResponseMetadata": {"HTTPStatusCode": 200}}

    def create_bucket(self, Bucket, **_):
        with self._lock:
            self._buckets.add(Bucket)
        return {"ResponseMetadata": {"HTTPStatusCode": 200}}

    def get_object(self, Bucket, Key, Range=None, **_):
        if self._is_data(Bucket, Key):
            if Range:
                rng = Range.replace("bytes=", "")
                left, right = rng.split("-")
                data = self._read_data_range(int(left), int(right))
            else:
                data = self._read_data_range(0, None)
            return {
                "Body": io.BytesIO(data),
                "ContentLength": len(data),
                "ResponseMetadata": {"HTTPStatusCode": 200},
            }
        with self._lock:
            store = self._meta.get(Bucket, {})
            if Key not in store:
                raise _client_error("NoSuchKey", "GetObject")
            data = store[Key]
        return {
            "Body": io.BytesIO(data),
            "ContentLength": len(data),
            "ResponseMetadata": {"HTTPStatusCode": 200},
        }

    def put_object(self, Bucket, Key, Body=None, **_):
        if hasattr(Body, "read"):
            data = Body.read()
        else:
            data = bytes(Body) if Body is not None else b""
        with self._lock:
            self._buckets.add(Bucket)
            self._meta.setdefault(Bucket, {})[Key] = data
        return {"ResponseMetadata": {"HTTPStatusCode": 200}}

    def delete_object(self, Bucket, Key, **_):
        with self._lock:
            self._meta.get(Bucket, {}).pop(Key, None)
        return {"ResponseMetadata": {"HTTPStatusCode": 200}}

    def upload_fileobj(self, Fileobj=None, Bucket=None, Key=None, **_):
        data = Fileobj.read()
        return self.put_object(Bucket=Bucket, Key=Key, Body=data)

    def upload_file(self, Filename=None, Bucket=None, Key=None, **_):
        with open(Filename, "rb") as f:
            data = f.read()
        return self.put_object(Bucket=Bucket, Key=Key, Body=data)


def make_cloud_object(sra_format_cls, fixture_path: str, accession: str,
                      bucket: str = "test", meta_bucket: str | None = None):
    """Build a CloudObject pointing at the local fixture using LocalFileS3.

    Returns the CloudObject. The boto3 PickleableS3ClientProxy is created
    inside a mock_aws block (so STS bootstrap succeeds) and is then
    swapped out for our LocalFileS3 stub before any real call is made.
    """
    from moto import mock_aws
    from dataplug import CloudObject

    if meta_bucket is None:
        meta_bucket = bucket + ".meta"

    storage = LocalFileS3(fixture_path, bucket, accession)
    storage.create_bucket(Bucket=meta_bucket)

    with mock_aws():
        co = CloudObject.from_s3(
            sra_format_cls,
            f"s3://{bucket}/{accession}",
            metadata_bucket=meta_bucket,
            fetch=False,
        )
    co._s3 = storage
    co.fetch()
    return co
