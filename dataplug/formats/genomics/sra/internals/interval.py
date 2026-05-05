from collections import defaultdict

Interval = tuple[int, int]


def ith_part(i: int, size: int, at: int = 0) -> Interval:
    return (at + size * i, at + size * (i + 1))


def partition_dry(total: int, num_chunks: int) -> list[Interval]:
    """ Divide range [0, total) into chunks of approximately equal size represented as
        tuples of [start, end).
    """
    if total == 0:
        return []

    num_chunks = max(1, num_chunks)
    num_chunks = min(total, num_chunks)

    chunk_size, leftovers = divmod(total, num_chunks)

    firsts = [ith_part(i, chunk_size + 1) for i in range(leftovers)]
    last_interval_end = firsts[-1][-1] if firsts else 0
    seconds = [ith_part(i, chunk_size, at=last_interval_end) for i in range(num_chunks - leftovers)]

    return firsts + seconds


def merge_intervals(intervals) -> list[Interval]:
    """ Given a list of intervals, merge overlapping intervals and return the merged list. """
    heights = defaultdict(lambda: 0)
    for tup in intervals:
        heights[tup[0]] += 1
        heights[tup[1]] -= 1

    merged = []
    height = 0
    left = min(heights.keys())
    for end in sorted(heights.keys()):
        height += heights[end]
        if left is None:
            left = end
        if height == 0:
            merged.append((left, end))
            left = None
    return merged
