import pytest

from dataplug.formats.genomics.sra.experimental import _walk_sample_tasks


def test_walk_sample_tasks_splits_samples_evenly():
    assert _walk_sample_tasks(sample_count=10, workers=3, overlap_steps=0) == [
        (0, 4),
        (4, 8),
        (8, 10),
    ]


def test_walk_sample_tasks_clamps_workers_to_sample_count():
    assert _walk_sample_tasks(sample_count=2, workers=8, overlap_steps=0) == [
        (0, 1),
        (1, 2),
    ]


def test_walk_sample_tasks_adds_clamped_overlap():
    assert _walk_sample_tasks(sample_count=10, workers=3, overlap_steps=2) == [
        (0, 6),
        (2, 10),
        (6, 10),
    ]


def test_walk_sample_tasks_handles_empty_input():
    assert _walk_sample_tasks(sample_count=0, workers=4, overlap_steps=0) == []


@pytest.mark.parametrize(
    ("sample_count", "workers", "overlap_steps"),
    [(-1, 1, 0), (1, 0, 0), (1, 1, -1)],
)
def test_walk_sample_tasks_rejects_invalid_values(sample_count, workers, overlap_steps):
    with pytest.raises(ValueError):
        _walk_sample_tasks(sample_count, workers, overlap_steps)
