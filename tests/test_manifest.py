from io import BytesIO
from unittest import TestCase

from mex_gene_archive.manifest import (
    compute_stream_md5sum,
)


class TestManifest(TestCase):
    def test_compute_stream_md5sum(self):
        zero_file = BytesIO(b"\x00" * 10000)
        digest = compute_stream_md5sum(zero_file).hexdigest()
        self.assertEqual("b85d6fb9ef4260dcf1ce0a1b0bff80d3", digest)
