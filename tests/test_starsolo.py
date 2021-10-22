import contextlib
import hashlib
from pathlib import Path
from unittest import TestCase

from mex_gene_archive.filter import (
    filter_mtx,
)
from mex_gene_archive.manifest import (
    compute_md5sums,
)
from mex_gene_archive.reader import (
    read_mex_archive,
    read_barcodes,
    read_features,
)
from mex_gene_archive.starsolo import (
    MULTIREAD_NAME,
    make_list_of_archive_files,
    make_output_type_term,
    validate_star_solo_out_arguments,
    archive_star_solo,
)
from .stage_data import srx5908538


def make_unique_filtered_archive(solo_dir, output_dir):
    config = {
        "experiment_accession": "SRP199641",
        "description": "scRNA-seq of HCC1395: LLU 10X Sequence1",
        "library_accession": "SRX5908538",
        "analysis_version": "version",
    }
    tar_name = archive_star_solo(
        solo_dir,
        config,
        "GeneFull_Ex50pAS",
        "Unique",
        "filtered",
        destination=output_dir,
    )
    return tar_name


class TestStarSolo(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.context = contextlib.ExitStack()
        cls.temp_dir = cls.context.enter_context(srx5908538())
        cls.analysis_dir = cls.temp_dir / "SRX5908538"
        cls.solo_dir = cls.analysis_dir / "Solo.out"
        cls.tar_name = make_unique_filtered_archive(cls.solo_dir, cls.analysis_dir)

    @classmethod
    def tearDownClass(cls):
        cls.context.pop_all()

    def test_make_list_of_archive_files_valid(self):
        root = Path("Solo.out")

        data = [
            ("Gene", "Unique", "raw"),
            ("Gene", "Unique", "filtered"),
            ("Gene", "EM", "raw"),
            ("Gene", "EM", "filtered"),
            ("GeneFull", "Unique", "raw"),
            ("GeneFull", "Unique", "filtered"),
            ("GeneFull", "EM", "raw"),
            ("GeneFull", "EM", "filtered"),
            ("GeneFull_Ex50pAS", "Unique", "raw"),
            ("GeneFull_Ex50pAS", "Unique", "filtered"),
            ("GeneFull_Ex50pAS", "EM", "raw"),
            ("GeneFull_Ex50pAS", "EM", "filtered"),
            ("SJ", "Unique", "raw"),
        ]

        for quantification, multiread, matrix in data:
            expected = make_list_of_archive_files(
                root, quantification, multiread, matrix
            )
            self.assertEqual(len(expected), 3)

            mex_matrix_parts = [
                "barcodes.tsv",
                "features.tsv",
                MULTIREAD_NAME[multiread],
            ]

            for i, pathname in enumerate(expected):
                self.assertEqual(pathname.parts[0], root.parts[0])
                self.assertEqual(pathname.parts[1], quantification)
                self.assertEqual(pathname.parts[2], matrix)
                self.assertEqual(pathname.name, mex_matrix_parts[i])

    def test_validate_star_solo_out_arguments_invalid(self):
        data = [
            ("gene", "Unique", "raw"),
            ("GeneFull", "EM", "filter"),
            ("GeneFull_Ex50pAS", "em", "filtered"),
            ("SJ", "EM", "raw"),
            ("SJ", "Unique", "filtered"),
        ]
        for quantification, multiread, matrix in data:
            self.assertRaises(
                ValueError,
                validate_star_solo_out_arguments,
                quantification,
                multiread,
                matrix,
            )

    def test_make_output_type_term(self):
        # Currently all gene quantification types have the same term
        expected_gene = {
            ("Unique", "filtered"): "sparse gene count matrix of unique reads",
            ("EM", "filtered"): "sparse gene count matrix of all reads",
            ("Unique", "raw"): "unfiltered sparse gene count matrix of unique reads",
            ("EM", "raw"): "unfiltered sparse gene count matrix of all reads",
        }
        gene_types = ["Gene", "GeneFull", "GeneFull_Ex50pAS"]
        for quantification in gene_types:
            for multiread, matrix in expected_gene:
                expected = expected_gene[(multiread, matrix)]
                term = make_output_type_term(quantification, multiread, matrix)
                self.assertEqual(expected, term)

        expected_splice = {
            (
                "SJ",
                "Unique",
                "raw",
            ): "unfiltered sparse splice junction count matrix of unique reads"
        }
        for quantification, multiread, matrix in expected_splice:
            expected = expected_splice[(quantification, multiread, matrix)]
            term = make_output_type_term(quantification, multiread, matrix)
            self.assertEqual(expected, term)

    def test_archive(self):
        expected = "GeneFull_Ex50pAS_Unique_filtered.tar.gz"
        self.assertTrue((self.analysis_dir / expected).is_file())
        filtered_dir = self.solo_dir / "GeneFull_Ex50pAS" / "filtered"

        data = read_mex_archive(self.tar_name)

        with open(filtered_dir / "matrix.mtx") as instream:
            for line in instream:
                if line.startswith("%"):
                    continue
                rows, columns, count = [int(x) for x in line.rstrip().split()]
                break

        self.assertEqual(data["matrix"].shape[0], rows)
        self.assertEqual(data["matrix"].shape[1], columns)

        with open(filtered_dir / "barcodes.tsv", "rt") as instream:
            barcodes = list(read_barcodes(instream))
            self.assertEqual(barcodes, data["barcodes"].to_list())

        with open(filtered_dir / "features.tsv", "rt") as instream:
            features = list(read_features(instream))
            gene_id = [x[0] for x in features]
            gene_symbols = [x[1] for x in features]
            self.assertEqual((len(features), 2), data["features"].shape)
            self.assertEqual(gene_id, data["features"]["gene_id"].to_list())
            self.assertEqual(gene_symbols, data["features"]["gene_symbols"].to_list())

    def test_filter(self):
        BARCODE_TSV_INDEX = 0
        FEATURE_TSV_INDEX = 1
        MTX_INDEX = 2

        raw_files = make_list_of_archive_files(
            self.solo_dir, "GeneFull_Ex50pAS", "Unique", "raw"
        )
        filtered_files = make_list_of_archive_files(
            self.solo_dir, "GeneFull_Ex50pAS", "Unique", "filtered"
        )

        files_md5 = compute_md5sums([filtered_files[MTX_INDEX]])

        md5 = hashlib.md5()
        for line in filter_mtx(
            raw_files[BARCODE_TSV_INDEX],
            raw_files[MTX_INDEX],
            filtered_files[BARCODE_TSV_INDEX],
        ):
            md5.update(line.encode("utf-8"))

        self.assertEqual(files_md5[0][1], md5.hexdigest())
