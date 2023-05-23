"""Functions for interacting with anndata matrices

Most likely
  write_anndata_as_mex_archive - write an AnnData matrix to a mex archive
  read_mex_archive_as_anndata - read a mex archive into an AnnData matrix
"""

from anndata import AnnData
import gzip
import hashlib
from io import BytesIO, StringIO
import os
from pathlib import Path
from scipy.io.mmio import MMFile
import tarfile

from mex_gene_archive.manifest import (
    ConfigError,
    create_metadata,
    write_manifest,
)
from mex_gene_archive.starsolo import (
    make_output_type_term,
    MULTIREAD_NAME,
    update_tarinfo,
    validate_star_solo_out_arguments,
)
from mex_gene_archive._version import __version__


class SimpleMMWriter(MMFile):
    # The default scipy writer defaults to more zeros
    # than we need.
    @staticmethod
    def _field_template(field, precision):
        return {
            MMFile.FIELD_REAL: "%.{}g\n".format(precision),
            MMFile.FIELD_INTEGER: "%i\n",
            MMFile.FIELD_UNSIGNED: "%u\n",
            MMFile.FIELD_COMPLEX: "%.{p}e %%.{p}e\n".format(p=precision),
        }.get(field, None)


def generate_barcode_bytes(adata):
    """ """
    for barcode in adata.obs_names:
        line = "{}{}".format(barcode, os.linesep)
        yield line.encode("ascii")


def generate_features_bytes(adata):
    for gene_id, gene_symbol in zip(adata.var_names, adata.var["gene_symbols"]):
        line = "{}\t{}\tGene Expression{}".format(gene_id, gene_symbol, os.linesep)
        yield line.encode("ascii")


def get_matrix_buffer(adata):
    with BytesIO() as stream:
        SimpleMMWriter().write(
            stream, adata.X.T, comment="", field=None, precision=None, symmetry=None
        )
        return stream.getvalue()


def compute_byte_md5sum(generator):
    md5 = hashlib.md5()
    for block in generator:
        md5.update(block)
    return md5.hexdigest()


def write_anndata_as_mex_archive(
    adata,
    quantification="GeneFull",
    multiread="Unique",
    matrix="raw",
    *,
    destination=None,
):
    """Write a mex archive from an anndata file.


    Parameters
    ----------
    """
    validate_star_solo_out_arguments(quantification, multiread, matrix)

    if "output_type" not in adata.uns:
        adata.uns["output_type"] = make_output_type_term(
            quantification, multiread, matrix
        )
    if "software_version" not in adata.uns:
        adata.uns["software_version"] = __version__
    if "arguments" not in adata.uns:
        adata.uns["arguments"] = ""

    required_terms = ["experiment_accession", "description", "library_accession"]
    for term in required_terms:
        if term not in adata.uns:
            raise ConfigError("required term {} missing from uns dictionary")

    archive_root = Path(quantification) / matrix

    manifest_filename = str(archive_root / "manifest.tsv")
    barcode_filename = str(archive_root / "barcodes.tsv")
    features_filename = str(archive_root / "features.tsv")
    matrix_filename = str(archive_root / MULTIREAD_NAME[multiread])

    tar_name = "{quantification}_{multiread}_{matrix}.tar.gz".format(
        quantification=quantification,
        multiread=multiread,
        matrix=matrix,
    )
    if destination is not None:
        tar_name = Path(destination) / tar_name

    md5s = [
        (barcode_filename, compute_byte_md5sum(generate_barcode_bytes(adata))),
        (features_filename, compute_byte_md5sum(generate_features_bytes(adata))),
        (matrix_filename, hashlib.md5(get_matrix_buffer(adata)).hexdigest()),
    ]
    manifest = create_metadata(adata.uns, md5s)

    manifest_buffer = BytesIO(
        write_manifest(StringIO(), manifest).getvalue().encode("ascii")
    )

    with gzip.GzipFile(tar_name, "wb", mtime=0) as gzipstream:
        with tarfile.open(
            mode="w", fileobj=gzipstream, format=tarfile.PAX_FORMAT
        ) as archive:
            info = tarfile.TarInfo(str(manifest_filename))
            update_tarinfo(info, fileobj=manifest_buffer)
            archive.addfile(info, manifest_buffer)

            with BytesIO() as stream:
                for line in generate_barcode_bytes(adata):
                    stream.write(line)
                stream.seek(0)
                info = tarfile.TarInfo(barcode_filename)
                update_tarinfo(info, fileobj=stream)
                archive.addfile(info, stream)

            with BytesIO() as stream:
                for line in generate_features_bytes(adata):
                    stream.write(line)
                stream.seek(0)
                info = tarfile.TarInfo(features_filename)
                update_tarinfo(info, fileobj=stream)
                archive.addfile(info, stream)

            with BytesIO() as stream:
                SimpleMMWriter().write(
                    stream,
                    adata.X.T,
                    comment="",
                    field=None,
                    precision=None,
                    symmetry=None,
                )
                stream.seek(0)

                info = tarfile.TarInfo(matrix_filename)
                update_tarinfo(info, fileobj=stream)
                archive.addfile(info, stream)

    return tar_name
