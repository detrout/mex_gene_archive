from io import StringIO
from itertools import combinations_with_replacement
import os
from unittest import TestCase

from mex_gene_archive.filter import (
    read_barcode_lineno_map,
    compute_raw_to_filtered_map,
)


def generate_fake_barcodes(length):
    bases = "AGTC"
    for barcode in combinations_with_replacement(bases, length):
        yield "".join(barcode)


class TestFilter(TestCase):
    # See also test_starsolo.TestStarSolo.test_filter for
    # a test the depends on available test data

    def test_read_barcode_lineno_map(self):
        fake_barcodes = StringIO(os.linesep.join(generate_fake_barcodes(4)))
        linenos = read_barcode_lineno_map(fake_barcodes)

        for i, barcode in enumerate(fake_barcodes):
            lineno = i + 1
            self.assertEqual(linenos[lineno], barcode)

    def test_compute_raw_to_filtered_map(self):
        """Generate a raw index to filtered index mapping"""

        filtered_fake = []
        for barcode in generate_fake_barcodes(4):
            if barcode[0] in ("A", "T"):
                filtered_fake.append(barcode)
        filtered_fake_barcodes = read_barcode_lineno_map(
            StringIO(os.linesep.join(filtered_fake))
        )
        raw_fake_barcodes = read_barcode_lineno_map(
            StringIO(os.linesep.join(generate_fake_barcodes(4)))
        )
        mapping = compute_raw_to_filtered_map(filtered_fake_barcodes, raw_fake_barcodes)
        print(mapping)
