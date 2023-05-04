import contextlib
import hashlib
from io import StringIO
import numpy
from pathlib import Path
from pytest import importorskip
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
    read_gene_features,
    read_sj_features,
    main as reader_main,
)
from mex_gene_archive.starsolo import (
    MULTIREAD_NAME,
    make_list_of_archive_files,
    make_output_type_term,
    make_tar_archive_name,
    validate_star_solo_out_arguments,
    parse_star_log_out_stream,
    archive_star_solo,
)
from .stage_data import (
    generate_count_matrix,
    make_sample_data,
    scratch_dir
)


def make_unique_filtered_archive(solo_dir, output_dir):
    config = {
        "experiment_accession": "SRP199641",
        "description": "scRNA-seq of HCC1395: LLU 10X Sequence1",
        "library_accession": "SRX5908538",
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
        cls.temp_dir = cls.context.enter_context(scratch_dir())
        print(cls.temp_dir)
        cls.analysis_dir = make_sample_data(cls.temp_dir)
        cls.solo_dir = cls.analysis_dir / "Solo.out"
        cls.tar_name = make_unique_filtered_archive(cls.solo_dir, cls.analysis_dir)

    @classmethod
    def tearDownClass(cls):
        cls.context.pop_all()

    def test_read_gene_features_null(self):
        instream = StringIO("100\tNULL\tGene Expression\n")
        features = read_gene_features(instream)
        self.assertEqual(features["gene_symbols"].iloc[0], "100")

    def test_read_sj_features(self):
        data = "chr1\t4878133\t4886743\t1\t1\t1\t0\t1\t30\n"
        instream = StringIO(data)
        expected = data.rstrip().split("\t")
        features = read_sj_features(instream)

        sj_columns = [
            "contig",
            "start",
            "end",
            "strand",
            "intron_motif",
            "annotated",
            "unique",
            "multi",
            "overhang",
        ]
        for i, name in enumerate(sj_columns):
            self.assertEqual(str(features[name].iloc[0]), expected[i])

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


    def test_parse_star_log_stream(self):
        software_version = "dev_EoI_2.7.9a_2021-09-30"
        arguments = """STAR --genomeDir genome --readFilesIn SRX5908538_R2.fastq.gz SRX5908538_R1.fastq.gz --readFilesCommand zcat --runThreadN 16 --genomeLoad NoSharedMemory --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMheaderCommentFile COfile.txt --outSAMheaderHD @HD VN:1.4 SO:coordinate --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD CB CR CY UB UR UY gx gn --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --sjdbScore 1 --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIdedup 1MM_CR --soloUMIfiltering MultiGeneUMI_CR --soloType CB_UMI_Simple --soloCellFilter EmptyDrops_CR --soloUMIlen 10 --soloCBlen 16 --soloBarcodeReadLength 0 --soloCBwhitelist 10xv2_allowlist.txt --soloStrand Forward --soloFeatures GeneFull_Ex50pAS SJ --soloMultiMappers Unique EM --limitBAMsortRAM 51539607552 --outTmpDir _STARtmp --outFileNamePrefix ./
"""

        head = StringIO("""STAR version={version}
STAR compilation time,server,dir=2021-09-30T17:07:55-07:00 :/tmp/tmp.01I63TZvRF/STAR/source
STAR git: On branch dev_ExonOverIntron ; commit 12beb0c52367d568fc993c1990795676e7f9d9cb ; diff files: 
##### Command Line:
{arguments}
##### Initial USER parameters from Command Line:
outFileNamePrefix                 ./
outTmpDir                         _STARtmp
###### All USER parameters from Command Line:
""".format(version=software_version, arguments=arguments))
        attributes = parse_star_log_out_stream(head)
        self.assertEqual(attributes["software_version"], software_version)
        self.assertEqual(attributes["arguments"], arguments)

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
            features = read_gene_features(instream)
            gene_id = features["gene_id"]
            gene_symbols = features["gene_symbols"]
            self.assertEqual((len(features), 2), data["features"].shape)
            self.assertEqual(gene_id.to_list(), data["features"]["gene_id"].to_list())
            self.assertEqual(gene_symbols.to_list(), data["features"]["gene_symbols"].to_list())

    def test_archive_to_h5ad(self):
        anndata = importorskip("anndata")
        target_file = self.temp_dir / "synthdata_GeneFull_Unique_filtered.h5ad"
        reader_main(["-o", str(target_file), str(self.tar_name)])

        self.assertTrue(target_file.is_file)

        adata = anndata.read_h5ad(target_file)

        filtered_dir = self.solo_dir / "GeneFull_Ex50pAS" / "filtered"
        with open(filtered_dir / "matrix.mtx") as instream:
            for line in instream:
                if line.startswith("%"):
                    continue
                rows, columns, count = [int(x) for x in line.rstrip().split()]
                break

        # AnnData matrix is transposed from what STAR writes
        self.assertEqual(adata.shape[0], columns)
        self.assertEqual(adata.shape[1], rows)

        with open(filtered_dir / "barcodes.tsv", "rt") as instream:
            barcodes = list(read_barcodes(instream))
            self.assertEqual(barcodes, adata.obs_names.to_list())

        with open(filtered_dir / "features.tsv", "rt") as instream:
            features = read_gene_features(instream)
            gene_id = features["gene_id"]
            gene_symbols = features["gene_symbols"]
            self.assertEqual(gene_id.to_list(), adata.var_names.to_list())
            self.assertEqual(gene_symbols.to_list(), adata.var["gene_symbols"].to_list())

    def test_generate_count_matrix(self):
        # Do we get the same count matrix for the same barcodes?
        count0 = generate_count_matrix(['AAAA', 'GGGG', 'TTTT'], 100)
        count1 = generate_count_matrix(['AAAA', 'TTTT'], 100)

        self.assertTrue(numpy.all(count0[:, 0] == count1[:, 0]))
        self.assertTrue(numpy.all(count0[:, 2] == count1[:, 1]))

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

    def test_make_tar_archive_name(self):
        dest = "/tmp"
        tar_name = make_tar_archive_name(self.solo_dir, "Gene", "Unique", "raw", dest)
        self.assertEqual(tar_name,  Path("/tmp/Gene_Unique_raw.tar.gz"))

        filename = "/tmp/archive.tar.gz"
        tar_name = make_tar_archive_name(
            self.solo_dir, "Gene", "Unique", "raw", filename)
        self.assertEqual(tar_name, Path(filename))

        filename = Path("/tmp/archive.tar.gz")
        tar_name = make_tar_archive_name(
            self.solo_dir, "Gene", "Unique", "raw", filename)
        self.assertEqual(tar_name, filename)

        tar_name = make_tar_archive_name(self.solo_dir, "Gene", "Unique", "raw", None)
        filename = self.solo_dir.parent / "Gene_Unique_raw.tar.gz"
        self.assertEqual(tar_name,  filename)

        self.assertRaises(
            RuntimeError, make_tar_archive_name, __file__, "Gene", "Unique", "raw", None)
