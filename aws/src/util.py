from urllib.parse import urlparse


def parse_s3_uri(uri: str) -> tuple[str, str]:
    """Parse an S3 URI into bucket and prefix.

    Args:
        uri (str): The S3 URI to parse.

    Returns:
        tuple[str, str]: A tuple containing the bucket name and the prefix.

    Raises:
        ValueError: If the URI scheme is not 's3'.
    """
    parsed = urlparse(uri)
    if parsed.scheme != "s3":
        raise ValueError(f"Invalid S3 URI: {uri}")
    bucket = parsed.netloc
    prefix = parsed.path.lstrip("/")
    return bucket, prefix
