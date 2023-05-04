from anndata import AnnData
from io import BytesIO
import numpy
from pathlib import Path
import pandas
from scipy.io import mmread
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

            adata.uns["experiment_accession"] = "test_experiment"
            adata.uns["description"] = "sample data generated for testing"
            adata.uns["library_accession"] = "test_library"

            write_anndata_as_mex_archive(
                adata, quantification, multiread, matrix, destination=tmp_dir
            )
            filename = tmp_dir / "{}_{}_{}.tar.gz".format(
                quantification, multiread, matrix
            )
            self.assertTrue(filename.exists())

            read = read_mex_archive_as_anndata(filename)

        self.assertTrue(numpy.all(adata.var_names == read.var_names))
        self.assertTrue(numpy.all(adata.obs_names == read.obs_names))
        self.assertTrue(numpy.all(adata.X == read.X))
