import pickle

import pytest

from dataplug.storage.picklableS3 import S3Path


def test_s3path_from_uri_exposes_bucket_key_and_uri():
    path = S3Path.from_uri("s3://example-bucket/foo/bar.txt")

    assert path.bucket == "example-bucket"
    assert path.key == "foo/bar.txt"
    assert path.as_uri() == "s3://example-bucket/foo/bar.txt"
    assert path.parts == ("/", "example-bucket", "foo", "bar.txt")


def test_s3path_from_bucket_key_strips_absolute_key():
    path = S3Path.from_bucket_key("example-bucket", "/foo/bar.txt")

    assert path.bucket == "example-bucket"
    assert path.key == "foo/bar.txt"


def test_s3path_rejects_relative_bucket_key_access():
    path = S3Path("foo/bar.txt")

    with pytest.raises(ValueError, match="relative path"):
        _ = path.bucket


def test_s3path_rejects_bucket_with_path_separator():
    with pytest.raises(ValueError, match="bucket argument contains"):
        S3Path.from_bucket_key("example-bucket/nested", "foo.txt")


def test_s3path_join_and_parent_reference_normalization():
    path = S3Path.from_bucket_key("example-bucket", "foo/old/../bar.txt")

    assert (path / "baz").key == "foo/bar.txt/baz"


def test_s3path_pickles_as_stable_value_object():
    path = S3Path.from_uri("s3://example-bucket/foo/bar.txt")

    assert pickle.loads(pickle.dumps(path)) == path


def test_s3path_behaves_like_posix_path_for_common_name_properties():
    path = S3Path.from_uri("s3://example-bucket/foo/bar.tar.gz")

    assert path.name == "bar.tar.gz"
    assert path.stem == "bar.tar"
    assert path.suffix == ".gz"
    assert path.suffixes == [".tar", ".gz"]


def test_s3path_parent_keeps_bucket_and_key_semantics():
    path = S3Path.from_uri("s3://example-bucket/foo/bar/baz.txt")

    assert path.parent.bucket == "example-bucket"
    assert path.parent.key == "foo/bar"
    assert path.parent.as_uri() == "s3://example-bucket/foo/bar"
    assert path.parents[0] == path.parent
    assert path.parents[1].key == "foo"


def test_s3path_parent_at_bucket_root_has_empty_key():
    path = S3Path.from_uri("s3://example-bucket/foo.txt")

    assert path.parent.bucket == "example-bucket"
    assert path.parent.key == ""
    assert path.parent.as_uri() == "s3://example-bucket/"


def test_s3path_joinpath_matches_division_operator():
    path = S3Path.from_bucket_key("example-bucket", "foo")

    assert path.joinpath("bar", "baz.txt") == path / "bar" / "baz.txt"
    assert path.joinpath("bar", "baz.txt").key == "foo/bar/baz.txt"


def test_s3path_relative_to_and_is_relative_to_accept_s3path_or_string():
    path = S3Path.from_uri("s3://example-bucket/foo/bar/baz.txt")

    assert path.relative_to(S3Path("/example-bucket/foo")).as_posix() == "bar/baz.txt"
    assert path.relative_to("/example-bucket/foo").as_posix() == "bar/baz.txt"
    assert path.is_relative_to(S3Path("/example-bucket/foo"))
    assert not path.is_relative_to("/other-bucket")


def test_s3path_with_name_stem_and_suffix_preserve_s3_type():
    path = S3Path.from_uri("s3://example-bucket/foo/bar.txt")

    assert path.with_name("baz.csv").as_uri() == "s3://example-bucket/foo/baz.csv"
    assert path.with_stem("baz").as_uri() == "s3://example-bucket/foo/baz.txt"
    assert path.with_suffix(".csv").as_uri() == "s3://example-bucket/foo/bar.csv"


def test_s3path_match_and_hash_behave_like_value_path():
    path = S3Path.from_uri("s3://example-bucket/foo/bar.txt")
    same_path = S3Path.from_bucket_key("example-bucket", "foo/bar.txt")

    assert path.match("*/foo/*.txt")
    assert {path: "value"}[same_path] == "value"
