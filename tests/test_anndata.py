from anndata import AnnData
from io import BytesIO, TextIOWrapper
import numpy
from pathlib import Path
import pandas
from scipy.io import mmread
import tarfile
import tempfile
from unittest import TestCase

from mex_gene_archive.anndata import (
    generate_barcode_bytes,
    generate_features_bytes,
    get_matrix_buffer,
    write_anndata_as_mex_archive,
)
from mex_gene_archive.manifest import ConfigError
from mex_gene_archive.reader import read_mex_archive_as_anndata
from mex_gene_archive.starsolo import MULTIREAD_NAME

from .stage_data import (
    generate_barcodes,
    generate_count_matrix,
    generate_features,
)


def make_adata(feature_length=20):
    barcodes = list(generate_barcodes())
    features = pandas.DataFrame(
        generate_features(20), columns=["id", "name", "type"]
    )
    count_matrix = generate_count_matrix(barcodes, feature_length).T
    adata = AnnData(X=count_matrix)
    adata.obs_names = barcodes
    adata.var_names = features["id"]
    adata.var["gene_symbols"] = features["name"]

    return adata


class TestMexAnndata(TestCase):
    def test_generate_barcode_bytes(self):
        adata = make_adata()

        expected_barcodes = "".join("{}\n".format(x) for x in adata.obs_names).encode(
            "ascii"
        )
        self.assertEqual(
            expected_barcodes, b"".join(x for x in generate_barcode_bytes(adata))
        )

    def test_generate_features_bytes(self):
        adata = make_adata()

        expected_features = []
        for gene_id, gene_name in zip(adata.var_names, adata.var["gene_symbols"]):
            expected_features.append(
                "\t".join((str(gene_id), str(gene_name), "Gene Expression\n"))
            )

        expected_features = "".join(expected_features).encode("ascii")
        self.assertEqual(
            expected_features, b"".join(x for x in generate_features_bytes(adata))
        )

    def test_get_matrix_buffer(self):
        adata = make_adata()

        matrix = mmread(BytesIO(get_matrix_buffer(adata)))

        # matrix market and anndata layouts are reversed
        self.assertTrue(numpy.all(adata.X.T == matrix))

    def test_write_and_read(self):
        adata = make_adata()

        quantification = "Gene"
        multiread = "Unique"
        matrix = "raw"

        expected_manifest = {
            "experiment_accession": "test_experiment",
            "description": "sample data generated for testing",
            "library_accession": "test_library",
        }

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_dir = Path(tmp_dir)
            self.assertRaises(
                ConfigError,
                write_anndata_as_mex_archive,
                adata,
                quantification,
                multiread,
                matrix,
                destination=tmp_dir,
            )

            for key in expected_manifest:
                adata.uns[key] = expected_manifest[key]

            write_anndata_as_mex_archive(
                adata, quantification, multiread, matrix, destination=tmp_dir
            )
            filename = tmp_dir / "{}_{}_{}.tar.gz".format(
                quantification, multiread, matrix
            )
            self.assertTrue(filename.exists())

            read = read_mex_archive_as_anndata(filename)

            with tarfile.open(name=filename, mode="r:*") as archive:
                for member in archive:
                    member_path = Path(member.name)
                    text_stream = TextIOWrapper(archive.extractfile(member), encoding="utf-8")
                    if member_path.name == "manifest.tsv":
                        metadata = pandas.read_csv(
                            text_stream, sep="\t", index_col="name")
                        break

        self.assertTrue(numpy.all(adata.var_names == read.var_names))
        self.assertTrue(numpy.all(adata.obs_names == read.obs_names))
        self.assertTrue(numpy.all(adata.X == read.X))

        for key in expected_manifest:
            self.assertEqual(adata.uns[key], expected_manifest[key])

        base = "{}/{}/".format(
            quantification, matrix
        )
        expected_md5s = {
            base + "barcodes.tsv": "md5sum:5c34cff881503b2434f6119bc223f3f6",
            base + "features.tsv": "md5sum:6f186211712a1568b0a58e44c9696ff3",
            base + MULTIREAD_NAME[multiread]: "md5sum:6dd76439cba2587e74ff7ba0b814fb6c",
        }
        for name in expected_md5s:
            self.assertEqual(metadata.loc[name, "value"], expected_md5s[name])
