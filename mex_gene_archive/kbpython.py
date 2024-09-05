import gzip
from io import BytesIO, StringIO
import json
import os
from pathlib import Path
import stat
import tarfile

from .manifest import (
    create_metadata, compute_md5sums, write_manifest)
from .utils import cast_to_list

KB_QUANTIFICATIONS = ["nascent", "mature", "ambiguous", "cell", "nucleus", "total"]
KB_RESULT_TYPES = ["counts_unfiltered", "counts_unfiltered_modified", "counts_filtered", "counts_filtered_modified"]


def archive_kbpython(root, config, quantifications, matrix, *, destination=None):
    """Archive a kbpython result directory into a tar file

    :param root: path to the directory that contains the result files
    :type root: str | Path
    :param config: metadata about the run, it really should include the source of the
                   raw data being processed, using values like experiment_accession for
                   ENCODE, or input_file_set for IGVF.
    :type config: dict
    :param quantifications: a name, or list of names for the matrices written by kbpython
                            ese KB_QUANTIFICATION for the known list
    :type quantifications: str | List[str]
    :param matrix: The name of the count directory used by kallisto. See KB_RESULT_TYPES
                   for the list of known values
    :type matrix: str
    :param destination: a filename or directory to write the output to. None will be treated
                        as the current directory. If destination is a directory, a filename
                        for the tar file will be computed based on the quantifications and
                        matrix parameters. If destination is a path to a file it wil be
                        used directly.
    :type destination: str | Path

    :returns: path to the output file
    :rtype: Path
    """
    quantifications = cast_to_list(quantifications)

    validate_kbpython_arguments(quantifications, matrix)

    archive_root = make_kbpython_root_name(root, quantifications, matrix)
    archive_files = make_list_of_archive_files(archive_root, quantifications, matrix)

    config["software"] = "kallisto"
    config.update(parse_kallisto_run_info(root / "run_info.json"))
    md5s = compute_md5sums(archive_files)
    manifest = create_metadata(config, md5s, root)

    manifest_buffer = BytesIO(
        write_manifest(StringIO(newline=""), manifest).getvalue().encode("utf-8")
    )
    manifest_filename = (archive_root / "manifest.tsv").relative_to(root)

    tar_name = make_tar_archive_name(root, quantifications, matrix, destination=destination)
    with gzip.GzipFile(tar_name, "wb", mtime=0) as gzipstream:
        with tarfile.open(mode="w", fileobj=gzipstream, format=tarfile.PAX_FORMAT) as archive:
            info = tarfile.TarInfo(str(manifest_filename))
            update_tarinfo(info, fileobj=manifest_buffer)
            archive.addfile(info, manifest_buffer)
            for filename in archive_files:
                info = tarfile.TarInfo(str(filename.relative_to(root)))
                update_tarinfo(info, filename)
                with open(filename, "rb") as instream:
                    archive.addfile(info, instream)

    return tar_name


def parse_kallisto_run_info(filename):
    """Open the kallisto run info file and read values

    :param filename: the path to a run_info.json file

    :returns: a dictionary of values. see
              parse_kallisto_run_info_stream for the details.

    """
    with open(filename, "rt") as instream:
        return parse_kallisto_run_info_stream(instream)

def parse_kallisto_run_info_stream(fileobj):
    """Open the kallisto run info file and read values

    :param fileobj: An input stream of a run_info.json file

    :returns: a dictionary of containing software, software_version,
              and arguments. read from the run_info file.
    """
    attributes = {}
    data = json.load(fileobj)
    attributes["software"] = "kallisto"
    attributes["software_version"] = data["kallisto_version"]
    attributes["arguments"] = data["call"]
    return attributes

def validate_kbpython_arguments(quantifications, matrix):
    """Make sure the arguments match the STAR Solo command line arguments

    Only certain combinations of arguments are possible.

    :param quantification: nascent, mature, ambiguous, cell, nucleus, total

    Raises
    ------

    ValueError if an invalid argument is provided
    """
    quantifications = cast_to_list(quantifications)

    for term in quantifications:
        if term not in KB_QUANTIFICATIONS:
            raise ValueError("Quantification type {} not in {}".format(quantification, KB_QUANTIFICATIONS))


    if matrix not in KB_RESULT_TYPES:
        raise ValueError("Matrix term {} not in {}".format(matrix, KB_RESULT_TYPES))

def make_quantification_name(quantifications):
    """Convert a list of quantifications into a name
    """
    print(quantifications)
    if len(quantifications) == 1 and quantifications[0] in KB_QUANTIFICATIONS:
        return quantifications[0]

    quantifications = set(quantifications)
    if quantifications == set(["nascent", "mature", "ambiguous"]):
        return "core"

    elif quantifications == set(["nascent", "mature", "ambiguous", "total"]):
        return "total+core"

    elif quantifications == set(KB_QUANTIFICATIONS):
        return "all"


    raise ValueError("Make up a name for {}".format(quantifications))


def make_kbpython_root_name(root, quantifications, matrix):
    """Compute the result directory from quantifications and matrix

    :param root: root directory for a kb_python run
    :param quantifications: list of quantification matrix types
    :param matrix: name of the type of count matrices.

    :returns: path to analysis directory for a specific matrix type
    """
    root = Path(root)
    validate_kbpython_arguments(quantifications, matrix)

    return root / matrix

def make_tar_archive_name(root, quantifications, matrix, destination):
    """Compute name of mex gene tar archives for different matrices
    """
    root = Path(root)
    quantifications = cast_to_list(quantifications)

    quantification_name = make_quantification_name(quantifications)
    tar_name = "{}_{}.tar.gz".format(matrix, quantification_name)
    if destination is not None:
        destination = Path(destination)
        if destination.is_dir():
            tar_name = destination / tar_name
        else:
            tar_name = destination
    elif root.is_dir():
        tar_name = root.parent / tar_name
    else:
        raise RuntimeError(
            "Unable to determine destination. Check {}, {}".format(destination, root))

    return tar_name

def make_list_of_archive_files(root, quantifications, matrix):
    """Generate list of files we expect to find in a UCI parse kb directory tree
    """
    quantifications = cast_to_list(quantifications)

    validate_kbpython_arguments(quantifications, matrix)
    root = Path(root)

    archive_root = make_kbpython_root_name(root, quantifications, matrix)

    archive_files = []
    for label in ["cells_x_genes.barcodes.txt", "cells_x_genes.genes.txt", "cells_x_genes.genes.names.txt"]:
        archive_files.append(root / label)
    for q in quantifications:
        archive_files.append(root / "cells_x_genes.{quantification}.mtx".format(quantification=q))

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
