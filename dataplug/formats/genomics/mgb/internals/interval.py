from __future__ import annotations

Interval = tuple[int, int]


def partition_dry(total: int, num_chunks: int) -> list[Interval]:
    if total == 0:
        return []

    num_chunks = max(1, min(total, num_chunks))
    chunk_size, leftovers = divmod(total, num_chunks)

    ranges = []
    start = 0
    for idx in range(num_chunks):
        size = chunk_size + (1 if idx < leftovers else 0)
        ranges.append((start, start + size))
        start += size
    return ranges


def merge_intervals(intervals) -> list[Interval]:
    sorted_intervals = sorted(intervals)
    if not sorted_intervals:
        return []

    merged = [sorted_intervals[0]]
    for left, right in sorted_intervals[1:]:
        prev_left, prev_right = merged[-1]
        if left <= prev_right:
            merged[-1] = (prev_left, max(prev_right, right))
        else:
            merged.append((left, right))
    return merged
