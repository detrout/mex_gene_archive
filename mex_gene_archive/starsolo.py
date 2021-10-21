from io import BytesIO, StringIO
import os
from pathlib import Path
import stat
import tarfile
import time

from .manifest import (
    compute_md5sums,
    create_metadata,
    write_manifest,
)


####
# functions for making archive file
MULTIREAD_NAME = {
    "Unique": "matrix.mtx",
    "Rescue": "UniqueAndMult-Rescue.mtx",
    "EM": "UniqueAndMult-EM.mtx",
}


def validate_star_solo_out_arguments(
    quantification="GeneFull", multiread="Unique", matrix="raw"
):
    quantification_terms = ["Gene", "GeneFull", "GeneFull_Ex50pAS", "SJ"]
    if quantification not in quantification_terms:
        raise ValueError("{} not in {}".format(quantification, quantification_terms))

    multiread_terms = ["Unique", "EM"]
    if multiread not in multiread_terms:
        raise ValueError("{} not in {}".format(multiread, multiread_terms))

    matrix_terms = ["filtered", "raw"]
    if matrix not in matrix_terms:
        raise ValueError("{} not in {}".format(matrix, matrix_terms))

    if quantification == "SJ":
        if multiread != "Unique":
            raise ValueError("Splice junctions do not support multread assignment")
        if matrix != "raw":
            raise ValueError("Splice junctions are only available as raw")


def make_list_of_archive_files(
    solo_root, quantification="GeneFull", multiread="Unique", matrix="raw"
):
    validate_star_solo_out_arguments(quantification, multiread, matrix)
    archive_files = []

    archive_files.append(solo_root / quantification / matrix / "barcodes.tsv")
    archive_files.append(solo_root / quantification / matrix / "features.tsv")

    archive_files.append(
        solo_root / quantification / matrix / MULTIREAD_NAME[multiread]
    )
    return archive_files


def update_tarinfo(info, filename):
    stat_info = os.stat(filename)
    info.size = stat_info[stat.ST_SIZE]
    info.mode = stat_info[stat.ST_MODE]
    info.mtime = time.time()
    info.uid = stat_info[stat.ST_UID]
    info.gid = stat_info[stat.ST_GID]
    info.type = tarfile.REGTYPE


def make_output_type_term(quantification="GeneFull", multiread="Unique", matrix="raw"):
    validate_star_solo_out_arguments(quantification, multiread, matrix)

    gene_term = {
        "Gene": "gene count matrix",
        "GeneFull": "gene count matrix",
        "GeneFull_Ex50pAS": "gene count matrix",
        "SJ": "splice junction count matrix",
    }[quantification]

    multiread_term = {
        "Unique": "unique",
        "EM": "all",
    }[multiread]

    matrix_term = {
        "filtered": "",
        "raw": "unfiltered ",
    }[matrix]

    output_type = "{count_matrix}sparse {quantification} of {multiread} reads".format(
        multiread=multiread_term,
        quantification=gene_term,
        count_matrix=matrix_term,
    )
    return output_type


def archive_star_solo(
    solo_root,
    config,
    quantification="GeneFull",
    multiread="Unique",
    matrix="raw",
    *,
    destination=None,
):
    validate_star_solo_out_arguments(quantification, multiread, matrix)

    archive_files = make_list_of_archive_files(
        solo_root, quantification, multiread, matrix
    )

    config['output_type'] = make_output_type_term(quantification, multiread, matrix)
    md5s = compute_md5sums(archive_files)
    manifest = create_metadata(config, md5s)
    manifest_buffer = BytesIO(
        write_manifest(StringIO(), manifest).getvalue().encode("utf-8")
    )

    tar_name = "{}_{}_{}.tar.gz".format(quantification, multiread, matrix)
    if destination is not None:
        tar_name = Path(destination) / tar_name
    elif solo_root.is_dir():
        tar_name = solo_root.parent / tar_name

    with tarfile.open(tar_name, "w:gz") as archive:
        info = tarfile.TarInfo("manifest.tsv")
        update_tarinfo(info, archive_files[0])
        info.size = len(manifest_buffer.getvalue())
        archive.addfile(info, manifest_buffer)
        for filename in archive_files:
            info = tarfile.TarInfo(str(filename.relative_to(solo_root)))
            update_tarinfo(info, filename)
            with open(filename, "rb") as instream:
                archive.addfile(info, instream)
