from __future__ import annotations

import logging
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Callable

    from dataplug.cloudobject import CloudObject


CHUNK_SIZE = 1048576
logger = logging.getLogger(__name__)


class DownloadError(Exception):
    pass


def get_status(resp) -> int:
    return resp.get("ResponseMetadata", {}).get("HTTPStatusCode")


def status_ok(response) -> bool:
    status = get_status(response)
    if status is None:
        return False
    status_success_low = 200
    status_success_high = 299
    return status_success_low <= status <= status_success_high


def stream_to_callback(stream, append: Callable[[bytes, int], None]) -> None:
    chunk = stream.read(CHUNK_SIZE)
    idx = 0
    while chunk != b"":
        append(chunk, idx)
        idx += len(chunk)
        chunk = stream.read(CHUNK_SIZE)
        logger.info("Next bytes %s", len(chunk))
    if hasattr(stream, "close"):
        stream.close()


def stream_to_memoryview(stream, destination: memoryview[int]) -> None:

    def memoryview_append(chunk: bytes, idx: int) -> None:
        destination[idx : idx + len(chunk)] = chunk

    stream_to_callback(stream, memoryview_append)


def stream_to_bytes(stream) -> bytes:
    chunks = []

    def bytes_append(chunk: bytes, idx: int) -> None:
        chunks.append(chunk)

    stream_to_callback(stream, bytes_append)
    return b"".join(chunks)


def download(cloud_object: CloudObject, destination: memoryview[int]) -> None:
    bucket = cloud_object.path.bucket
    key = cloud_object.path.key
    storage = cloud_object.storage
    resp = storage.get_object(Bucket=bucket, Key=key)
    if not status_ok(resp):
        raise DownloadError
    data_stream = resp["Body"]

    stream_to_memoryview(data_stream, destination)


def download_range(cloud_object, left, right):
    bucket = cloud_object.path.bucket
    key = cloud_object.path.key
    res = cloud_object.storage.get_object(
        Bucket=bucket,
        Key=key,
        Range=f"bytes={left}-{right}",
    )
    return stream_to_bytes(res["Body"])
