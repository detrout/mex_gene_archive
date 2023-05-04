import gzip
from io import BytesIO, StringIO
import os
from pathlib import Path
import stat
import tarfile

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
    """Make sure the arguments match the STAR Solo command line arguments
    """
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


def make_archive_root_name(solo_root, quantification, multiread, matrix):
    """Generate list of files we expect to find in a STAR Solo.out directory tree
    """
    validate_star_solo_out_arguments(quantification, multiread, matrix)

    return solo_root / quantification / matrix


def make_tar_archive_name(
        solo_root, quantification, multiread, matrix, destination):
    """Compute name of mex gene tar archives for different matrices
    """
    solo_root = Path(solo_root)
    tar_name = "{}_{}_{}.tar.gz".format(quantification, multiread, matrix)
    if destination is not None:
        destination = Path(destination)
        if destination.is_dir():
            tar_name = destination / tar_name
        else:
            tar_name = destination
    elif solo_root.is_dir():
        tar_name = solo_root.parent / tar_name
    else:
        raise RuntimeError(
            "Unable to determine destination. Check {}, {}".format(destination, solo_root))

    return tar_name


def make_list_of_archive_files(
    solo_root, quantification="GeneFull", multiread="Unique", matrix="raw"
):
    """Generate list of files we expect to find in a STAR Solo.out directory tree
    """
    validate_star_solo_out_arguments(quantification, multiread, matrix)
    archive_files = []

    archive_root = make_archive_root_name(
        solo_root, quantification, multiread, matrix)

    archive_files.append(archive_root / "barcodes.tsv")
    archive_files.append(archive_root / "features.tsv")

    archive_files.append(
        archive_root / MULTIREAD_NAME[multiread]
    )
    return archive_files


def update_tarinfo(info, filename=None, fileobj=None):
    """Fill in tarinfo fields for making tar archive
    """
    if filename is not None:
        stat_info = os.stat(filename)
        info.size = stat_info[stat.ST_SIZE]
    elif fileobj is not None:
        curpos = fileobj.tell()
        size = fileobj.seek(0, 2)
        fileobj.seek(curpos, 0)
        info.size = size

    info.mtime = 0
    info.mode = 0o644
    info.uid = 0
    info.gid = 0
    info.uname = "root"
    info.gname = "root"
    info.type = tarfile.REGTYPE


def make_output_type_term(quantification="GeneFull", multiread="Unique", matrix="raw"):
    """Generate ENCODE controlled vocabulary term for the specified result type

    The different combinations of quantification, multiread, and matrix are represented
    as different object types and the ENCODE portal.
    """
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


def parse_star_log_out(filename):
    """Wrapper for parse_star_log_out_stream to do file IO
    """
    with open(filename, "rt") as instream:
        return parse_star_log_out_stream(instream)


def parse_star_log_out_stream(fileobj):
    """Read a stream holding STAR Log.out and return version and arguments
    """
    star_version_prefix = "STAR version="
    attributes = {}
    for line in fileobj:
        if line.startswith(star_version_prefix):
            attributes["software_version"] = line.rstrip()[len(star_version_prefix):]
        elif line.startswith("##### Command Line:"):
            attributes["arguments"] = next(fileobj)

    return attributes


def archive_star_solo(
    solo_root,
    config,
    quantification="GeneFull",
    multiread="Unique",
    matrix="raw",
    *,
    destination=None,
):
    """Archive a specific STAR solo result directory


    Parameters
    ----------
    solo_root : path to STAR's Solo.out directory where
        a file named {quantification}_{multiread}_{matrix}.tar.gz will be
        written.
    config : dictionary of configuration options
    quantification : Which counting method to use "Gene", "GeneFull",
        "GeneFull_Ex50pAS"
    multiread : which STAR EM processing level to use "Unique", "EM"
    matrix : which matrix to read either "raw" or "filtered"
    destination : what directory to write the archive to, defaults to
        solo_root/..
    """
    validate_star_solo_out_arguments(quantification, multiread, matrix)

    archive_files = make_list_of_archive_files(
        solo_root, quantification, multiread, matrix
    )

    config['output_type'] = make_output_type_term(quantification, multiread, matrix)
    config.update(parse_star_log_out(solo_root / ".." / "Log.out"))
    md5s = compute_md5sums(archive_files)
    manifest = create_metadata(config, md5s)
    manifest_buffer = BytesIO(
        write_manifest(StringIO(), manifest).getvalue().encode("utf-8")
    )
    archive_root = make_archive_root_name(
        solo_root, quantification, multiread, matrix)
    manifest_filename = (archive_root / "manifest.tsv").relative_to(solo_root)

    tar_name = make_tar_archive_name(
        solo_root, quantification, multiread, matrix, destination)
    with gzip.GzipFile(tar_name, "wb", mtime=0) as gzipstream:
        with tarfile.open(mode="w", fileobj=gzipstream, format=tarfile.PAX_FORMAT) as archive:
            info = tarfile.TarInfo(str(manifest_filename))
            update_tarinfo(info, fileobj=manifest_buffer)
            archive.addfile(info, manifest_buffer)
            for filename in archive_files:
                info = tarfile.TarInfo(str(filename.relative_to(solo_root)))
                update_tarinfo(info, filename)
                with open(filename, "rb") as instream:
                    archive.addfile(info, instream)

    return tar_name
