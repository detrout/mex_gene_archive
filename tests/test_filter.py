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

        # The generate fake barcodes never gets to C and only includes the following
        # G lines.
        # 'ACCC': 20,
        # 'GGGG': 21,
        # 'GGGT': 22,
        # 'GGGC': 23,
        # 'GGTT': 24,
        # 'GGTC': 25,
        # 'GGCC': 26,
        # 'GTTT': 27,
        # 'GTTC': 28,
        # 'GTCC': 29,
        # 'GCCC': 30,
        # 'TTTT': 31,
        # So it's lines 21-30 that get deleted and need to be remapped.

        mapping = compute_raw_to_filtered_map(filtered_fake_barcodes, raw_fake_barcodes)
        expected_mapping = {
            1: 1,
            2: 2,
            3: 3,
            4: 4,
            5: 5,
            6: 6,
            7: 7,
            8: 8,
            9: 9,
            10: 10,
            11: 11,
            12: 12,
            13: 13,
            14: 14,
            15: 15,
            16: 16,
            17: 17,
            18: 18,
            19: 19,
            20: 20,
            31: 21,
            32: 22,
            33: 23,
            34: 24,
        }

        self.assertEqual(mapping, expected_mapping)
