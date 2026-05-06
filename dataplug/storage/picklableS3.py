from __future__ import annotations

import json
import logging
import os
import time
import uuid
from contextlib import suppress
from pathlib import PurePosixPath
from typing import TYPE_CHECKING

import boto3
import botocore.client

if TYPE_CHECKING:

    from mypy_boto3_s3.type_defs import DeleteObjectOutputTypeDef, DeleteObjectsOutputTypeDef

logger = logging.getLogger(__name__)

S3_FULL_ACCESS_POLICY = json.dumps(
    {
        "Id": "BucketPolicy",
        "Version": "2012-10-17",
        "Statement": [
            {
                "Sid": "AllAccess",
                "Action": "s3:*",
                "Effect": "Allow",
                "Resource": ["arn:aws:s3:::*", "arn:aws:s3:::*/*"],
                "Principal": "*",
            }
        ],
    }
)


class PickleableS3ClientProxy:
    """
    A Pickleable S3 client proxy that can be pickled and unpickled.

    Dataplug requires having an S3 client that can be pickled and sent remotely to workers, so that
    remote workers can access S3 objects. This class is a proxy to an S3 client that can be pickled.
    For it, we request temporary credentials to S3 using STS. To request temporary credentials, we
    authenticate using the locally configured credentials (through the environment variables or the
    AWS configuration file) and assume a role with full access to S3. The temporary credentials are
    saved in the instance and used to create an S3 client, and to create a new client when the object
    is unpickled.
    """

    def __init__(
        self,
        region_name: str | None = None,
        endpoint_url: str | None = None,
        credentials: dict | None = None,
        role_arn: str | None = None,
        token_duration_seconds: int | None = None,
        botocore_config_kwargs: dict | None = None,
    ):
        self.region_name = region_name
        self.endpoint_url = endpoint_url
        self.botocore_config_kwargs = botocore_config_kwargs or {}
        self.role_arn = role_arn
        self.session_name = None
        self.token_duration_seconds = token_duration_seconds or 86400  # 24 hours

        logger.debug("Requesting temporary credentials for S3 authentication")
        sts_admin = boto3.client(
            "sts",
            region_name=region_name,
            endpoint_url=endpoint_url,
            config=botocore.client.Config(**self.botocore_config_kwargs),
        )

        if credentials:
            # check if credentials are valid
            if not all(key in credentials for key in ["AccessKeyId", "SecretAccessKey"]):
                raise ValueError(
                    "Invalid credentials. AccessKeyId and SecretAccessKey are required if credentials are provided."
                )
            self.credentials = credentials

        elif role_arn is not None:
            self.session_name = "-".join(["dataplug", str(int(time.time())), uuid.uuid4().hex])
            logger.debug(
                "Assuming role %s with generated session name %s",
                self.role_arn,
                self.session_name,
            )

            response = sts_admin.assume_role(
                RoleArn=self.role_arn,
                RoleSessionName=self.session_name,
                Policy=S3_FULL_ACCESS_POLICY,
                DurationSeconds=self.token_duration_seconds,
            )
            self.credentials = response["Credentials"]

        else:
            logger.debug("Getting session token")
            response = sts_admin.get_session_token(DurationSeconds=self.token_duration_seconds)
            self.credentials = response["Credentials"]

        self.__client = boto3.client(
            "s3",
            aws_access_key_id=self.credentials["AccessKeyId"],
            aws_secret_access_key=self.credentials["SecretAccessKey"],
            aws_session_token=self.credentials.get("SessionToken"),
            endpoint_url=self.endpoint_url,
            region_name=self.region_name,
            config=botocore.client.Config(**self.botocore_config_kwargs),
        )

    def _new_client(self):
        session = boto3.Session(
            aws_access_key_id=self.credentials["AccessKeyId"],
            aws_secret_access_key=self.credentials["SecretAccessKey"],
            aws_session_token=self.credentials.get("SessionToken"),
            region_name=self.region_name,
        )
        return session.client(
            "s3",
            endpoint_url=self.endpoint_url,
            config=botocore.client.Config(**self.botocore_config_kwargs),
        )

    def __getstate__(self):
        logger.debug("Pickling S3 client")
        return {
            "credentials": self.credentials,
            "endpoint_url": self.endpoint_url,
            "region_name": self.region_name,
            "botocore_config_kwargs": self.botocore_config_kwargs,
            "role_arn": self.role_arn,
            "session_name": self.session_name,
            "token_duration_seconds": self.token_duration_seconds,
        }

    def __setstate__(self, state):
        logger.debug("Restoring S3 client")
        self.credentials = state["credentials"]
        self.endpoint_url = state["endpoint_url"]
        self.region_name = state["region_name"]
        self.botocore_config_kwargs = state["botocore_config_kwargs"]
        self.role_arn = state["role_arn"]
        self.session_name = state["session_name"]
        self.token_duration_seconds = state["token_duration_seconds"]

        self.__client = boto3.client(
            "s3",
            aws_access_key_id=self.credentials["AccessKeyId"],
            aws_secret_access_key=self.credentials["SecretAccessKey"],
            aws_session_token=self.credentials.get("SessionToken"),
            endpoint_url=self.endpoint_url,
            region_name=self.region_name,
            config=botocore.client.Config(**self.botocore_config_kwargs),
        )

    def abort_multipart_upload(self, *args, **kwargs):
        response = self.__client.abort_multipart_upload(*args, **kwargs)
        logger.debug("%s", response.get("ResponseMetadata", {}))
        return response

    def complete_multipart_upload(self, *args, **kwargs):
        response = self.__client.complete_multipart_upload(*args, **kwargs)
        logger.debug("%s", response.get("ResponseMetadata", {}))
        return response

    def create_multipart_upload(self, *args, **kwargs):
        response = self.__client.create_multipart_upload(*args, **kwargs)
        logger.debug("%s", response.get("ResponseMetadata", {}))
        return response

    def download_file(self, *args, **kwargs):
        self.__client.download_file(*args, **kwargs)
        logger.debug("%s", {"HTTP-Status": "200"})

    def download_fileobj(self, *args, **kwargs):
        self.__client.download_fileobj(*args, **kwargs)
        logger.debug("%s", {"HTTP-Status": "200"})

    def generate_presigned_post(self, *args, **kwargs):
        response = self.__client.generate_presigned_post(*args, **kwargs)
        logger.debug("%s", response.get("ResponseMetadata", {}))
        return response

    def generate_presigned_url(self, *args, **kwargs):
        response = self.__client.generate_presigned_url(*args, **kwargs)
        logger.debug("%s", response)
        return response

    def get_object(self, *args, **kwargs):
        response = self.__client.get_object(*args, **kwargs)
        logger.debug("%s", response.get("ResponseMetadata", {}))
        return response

    def delete_object(self, *args, **kwargs) -> DeleteObjectOutputTypeDef:
        response = self.__client.delete_object(*args, **kwargs)
        logger.debug("%s", response.get("ResponseMetadata", {}))
        return response

    def delete_objects(self, *args, **kwargs) -> DeleteObjectsOutputTypeDef:
        response = self.__client.delete_objects(*args, **kwargs)
        logger.debug("%s", response.get("ResponseMetadata", {}))
        return response

    def head_bucket(self, *args, **kwargs):
        response = self.__client.head_bucket(*args, **kwargs)
        logger.debug("%s", response.get("ResponseMetadata", {}))
        return response

    def head_object(self, *args, **kwargs):
        response = self.__client.head_object(*args, **kwargs)
        logger.debug("%s", response.get("ResponseMetadata", {}))
        return response

    def list_buckets(self):
        response = self.__client.list_buckets()
        logger.debug("%s", response.get("ResponseMetadata", {}))
        return response

    def list_multipart_uploads(self, *args, **kwargs):
        response = self.__client.list_multipart_uploads(*args, **kwargs)
        logger.debug("%s", response.get("ResponseMetadata", {}))
        return response

    def list_objects(self, *args, **kwargs):
        response = self.__client.list_objects(*args, **kwargs)
        logger.debug("%s", response.get("ResponseMetadata", {}))
        return response

    def list_objects_v2(self, *args, **kwargs):
        response = self.__client.list_objects_v2(*args, **kwargs)
        logger.debug("%s", response.get("ResponseMetadata", {}))
        return response

    def list_parts(self, *args, **kwargs):
        response = self.__client.list_parts(*args, **kwargs)
        logger.debug("%s", response.get("ResponseMetadata", {}))
        return response

    def put_object(self, *args, **kwargs):
        response = self.__client.put_object(*args, **kwargs)
        logger.debug("%s", response.get("ResponseMetadata", {}))
        return response

    def upload_file(self, *args, **kwargs):
        self.__client.upload_file(*args, **kwargs)
        logger.debug("%s", {"HTTP-Status": "200"})

    def upload_fileobj(self, *args, **kwargs):
        self.__client.upload_fileobj(*args, **kwargs)
        logger.debug("%s", {"HTTP-Status": "200"})

    def upload_part(self, *args, **kwargs):
        response = self.__client.upload_part(*args, **kwargs)
        logger.debug("%s", response.get("ResponseMetadata", {}))
        return response

    def create_bucket(self, *args, **kwargs):
        response = self.__client.create_bucket(*args, **kwargs)
        logger.debug("%s", response.get("ResponseMetadata", {}))
        return response


class S3Path:
    """
    Small POSIX-like path value for AWS S3 object references.

    pathlib internals changed across Python versions, so this class wraps the
    public PurePosixPath API instead of subclassing PurePath or overriding
    private path flavours.
    """

    __slots__ = ("_path",)

    def __init__(self, *parts: object):
        if not parts:
            raw_path = "."
        else:
            raw_path = "/".join(str(part).strip("/") for part in parts if str(part))
            if str(parts[0]).startswith("/"):
                raw_path = "/" + raw_path
        self._path = PurePosixPath(self._normalize(raw_path))

    @staticmethod
    def _normalize(path: str) -> str:
        is_absolute = path.startswith("/")
        normalized_parts: list[str] = []
        for part in path.split("/"):
            if part in ("", "."):
                continue
            if part == "..":
                if normalized_parts and normalized_parts[-1] != "..":
                    normalized_parts.pop()
                elif not is_absolute:
                    normalized_parts.append(part)
                continue
            normalized_parts.append(part)

        normalized = "/".join(normalized_parts)
        if is_absolute:
            return "/" + normalized if normalized else "/"
        return normalized or "."

    @classmethod
    def from_uri(cls, uri: str) -> S3Path:
        """
        from_uri class method create a class instance from url

        >> from s3path import PureS3Path
        >> PureS3Path.from_url('s3://<bucket>/<key>')
        << PureS3Path('/<bucket>/<key>')
        """
        if not uri.startswith("s3://"):
            raise ValueError("Provided uri seems to be no S3 URI!")
        return cls("/", uri[5:])

    @classmethod
    def from_bucket_key(cls, bucket: str, key: str) -> S3Path:
        """
        from_bucket_key class method create a class instance from bucket, key pair's

        >> from s3path import PureS3Path
        >> PureS3Path.from_bucket_key(bucket='<bucket>', key='<key>')
        << PureS3Path('/<bucket>/<key>')
        """
        bucket = cls("/", bucket)
        if len(bucket.parts) != 2:
            raise ValueError(
                f"bucket argument contains more then one path element: {bucket}"
            )
        key = cls(key)
        if key.is_absolute():
            key = key.relative_to("/")
        return bucket / key

    @property
    def bucket(self) -> str:
        """
        The AWS S3 Bucket name, or ''
        """
        self._absolute_path_validation()
        with suppress(ValueError):
            _, bucket, *_ = self.parts
            return bucket
        return ""

    @property
    def key(self) -> str:
        """
        The AWS S3 Key name, or ''
        """
        self._absolute_path_validation()
        return "/".join(self.parts[2:])

    @property
    def virtual_directory(self) -> str:
        """
        The parent virtual directory of a key
        Example: foo/bar/baz -> foo/baz
        """
        vdir, _ = self.key.rsplit("/", 1)
        return vdir

    def as_uri(self) -> str:
        """
        Return the path as a 's3' URI.
        """
        self._absolute_path_validation()
        return f"s3://{self.bucket}/{self.key}"

    @property
    def parts(self) -> tuple[str, ...]:
        return self._path.parts

    @property
    def name(self) -> str:
        return self._path.name

    @property
    def suffix(self) -> str:
        return self._path.suffix

    @property
    def suffixes(self) -> list[str]:
        return self._path.suffixes

    @property
    def stem(self) -> str:
        return self._path.stem

    @property
    def parent(self) -> S3Path:
        return self.__class__(self._path.parent)

    @property
    def parents(self) -> tuple[S3Path, ...]:
        return tuple(self.__class__(parent) for parent in self._path.parents)

    def is_absolute(self) -> bool:
        return self._path.is_absolute()

    def relative_to(self, other: str | S3Path) -> S3Path:
        other_path = other._path if isinstance(other, S3Path) else PurePosixPath(str(other))
        return self.__class__(self._path.relative_to(other_path))

    def is_relative_to(self, other: str | S3Path) -> bool:
        other_path = other._path if isinstance(other, S3Path) else PurePosixPath(str(other))
        return self._path.is_relative_to(other_path)

    def as_posix(self) -> str:
        return self._path.as_posix()

    def joinpath(self, *other: object) -> S3Path:
        return self.__class__(self._path, *other)

    def match(self, path_pattern: str) -> bool:
        return self._path.match(path_pattern)

    def with_name(self, name: str) -> S3Path:
        return self.__class__(self._path.with_name(name))

    def with_stem(self, stem: str) -> S3Path:
        return self.__class__(self._path.with_stem(stem))

    def with_suffix(self, suffix: str) -> S3Path:
        return self.__class__(self._path.with_suffix(suffix))

    def __truediv__(self, other: object) -> S3Path:
        return self.joinpath(other)

    def __str__(self) -> str:
        return self.as_posix()

    def __fspath__(self) -> str:
        return self.as_posix()

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, S3Path):
            return NotImplemented
        return self._path == other._path

    def __hash__(self) -> int:
        return hash(self._path)

    def _absolute_path_validation(self):
        if not self.is_absolute():
            raise ValueError("relative path have no bucket, key specification")

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(bucket={self.bucket},key={self.key})"


os.environ.setdefault("AWS_REQUEST_CHECKSUM_CALCULATION", "WHEN_REQUIRED")
