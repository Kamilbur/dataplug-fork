from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Callable
    from dataplug.cloudobject import CloudObject


CHUNK_SIZE = 1048576


def stream_to_callback(stream, append: Callable[[bytes, int], None]) -> None:
    chunk = stream.read(CHUNK_SIZE)
    idx = 0
    while chunk != b"":
        append(chunk, idx)
        idx += len(chunk)
        chunk = stream.read(CHUNK_SIZE)
    if hasattr(stream, "close"):
        stream.close()


def stream_to_bytes(stream) -> bytes:
    chunks = []

    def append(chunk: bytes, idx: int) -> None:
        chunks.append(chunk)

    stream_to_callback(stream, append)
    return b"".join(chunks)


def download_range(cloud_object: CloudObject, left: int, right: int) -> bytes:
    response = cloud_object.storage.get_object(
        Bucket=cloud_object.path.bucket,
        Key=cloud_object.path.key,
        Range=f"bytes={left}-{right}",
    )
    return stream_to_bytes(response["Body"])
