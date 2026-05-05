from __future__ import annotations

import io
import logging
from pathlib import Path
from typing import TYPE_CHECKING

from dataplug.entities import CloudDataFormat, CloudObjectSlice, PartitioningStrategy
from dataplug.preprocessing.metadata import PreprocessingMetadata

if TYPE_CHECKING:
    from dataplug.cloudobject import CloudObject

    from .internals.sra import SRALines

logger = logging.getLogger(__name__)


def preprocess_sra(cloud_object: CloudObject):
    from .internals.sra import VColumns

    logger.info("Preprocessing sra started")

    acc = Path(cloud_object.path.as_uri()).name.split(".")[0]
    with VColumns.from_filepath(acc) as vcols:
        num_lines = len(vcols)
    return PreprocessingMetadata(
        metadata=io.BytesIO(b"nempty"),
        attributes={
            "total_lines": num_lines,
        },
    )


@CloudDataFormat(preprocessing_function=preprocess_sra)
class SRA:
    total_lines: int
    index_key: str


class SRASlice(CloudObjectSlice):
    def __init__(self, start, end, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.start = start
        self.end = end

    def get(self, split: bool = True) -> SRALines:
        return list(self.lazy_get(split=split))

    def lazy_get(self, split: bool = True, write_unpaired: bool = False):
        from .internals.sra import VColumns, _format_spot

        if not self.cloud_object:
            msg = "Cloud object is not set for this slice."
            raise ValueError(msg)
        acc = self.cloud_object.path.as_uri()
        with VColumns.from_filepath(acc) as vcols:
            for row_idx in range(self.start + 1, self.end + 1):
                row = [col.read(row_idx) for col in vcols.columns]
                yield _format_spot(acc, row_idx, *row, split=split)


def partition_into_ranges(total_lines: int, num_chunks: int) -> list[tuple[int, int]]:
    n, r = divmod(total_lines, num_chunks)

    if total_lines == 0:
        return []
    if total_lines < num_chunks:
        logger.warning(
            "Number of chunks given is greater "
            "than total number of lines in file. "
            "Dividing files into single-line chunks... "
            "For performance gains, consider reducing "
            "number of chunks."
        )
        num_chunks = total_lines

    n, r = divmod(total_lines, num_chunks)  # lines_per_chunk, leftovers
    firsts = [((n + 1) * i, (n + 1) * (i + 1)) for i in range(r)]
    last = firsts[-1][-1] if firsts else 0
    seconds = [(last + n * i, last + n * (i + 1)) for i in range(num_chunks - r)]
    return firsts + seconds


@PartitioningStrategy(dataformat=SRA)
def partition_chunks_strategy(cloud_object: CloudObject, num_chunks: int) -> list[SRASlice]:
    logger.info("SRA partitioning started")
    total_lines = int(cloud_object.get_attribute("total_lines"))

    ranges = partition_into_ranges(total_lines, num_chunks)

    return [SRASlice(*r) for r in ranges]
