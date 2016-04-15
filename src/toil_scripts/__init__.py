import logging

_log = logging.getLogger(__name__)

# TODO: replace this with toil_scripts.lib.urls._download_s3_url

def download_from_s3_url(file_path, url):
    from urlparse import urlparse
    from boto.s3.connection import S3Connection
    s3 = S3Connection()
    try:
        parsed_url = urlparse(url)
        if not parsed_url.netloc or not parsed_url.path.startswith('/'):
            raise RuntimeError("An S3 URL must be of the form s3:/BUCKET/ or "
                               "s3://BUCKET/KEY. '%s' is not." % url)
        bucket = s3.get_bucket(parsed_url.netloc)
        key = bucket.get_key(parsed_url.path[1:])
        _log.info("From URL %s, parsed --> %s, bucket %s, key %s", url, parsed_url, bucket, key)
        key.get_contents_to_filename(file_path)
    finally:
        s3.close()
