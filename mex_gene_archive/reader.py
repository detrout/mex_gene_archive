import argparse
import csv
from io import TextIOWrapper
import logging
from pathlib import Path
import pandas
from scipy.io import mmread
import tarfile

from .manifest import read_manifest


logger = logging.getLogger(__name__)


def read_barcodes(instream):
    for line in instream:
        yield line.rstrip()


def read_features(instream):
    reader = csv.reader(instream, delimiter="\t")
    for row in reader:
        # Third column is frequently boring, why do we even save it?
        yield (row[0], row[1])


def remove_manifest_md5s(manifest):
    to_delete = []
    for key in manifest:
        path_name = Path(key)
        if (
            path_name.name in ("manifest.tsv", "barcodes.tsv", "features.tsv")
            or path_name.suffix == ".mtx"
        ):
            to_delete.append(key)

    for key in to_delete:
        del manifest[key]

    return manifest


def read_mex_archive(filename=None, fileobj=None):
    """Reads a mex archive into a dictionary

    Parameters
    ----------
    filename : Name of file to read
    fileobj : file-like object to read

    one and only one of filename or fileobj needs to be provided.

    Returns
    -------
    A dictionary with the manifest, (cell) barcodes, (gene) features,
    and the sparse matrix in features in rows and and barcodes in columns
    format.

    note that the preferred AnnData format is transposed from the default
    mapper format.
    """
    result = {}

    with tarfile.open(name=filename, mode="r:*", fileobj=fileobj) as archive:
        for member in archive:
            member_path = Path(member.name)
            text_stream = TextIOWrapper(archive.extractfile(member), encoding="utf-8")
            if member_path.name == "manifest.tsv":
                result["metadata"] = read_manifest(text_stream)
            elif member_path.name == "barcodes.tsv":
                result["barcodes"] = pandas.Series(read_barcodes(text_stream))
            elif member_path.name == "features.tsv":
                result["features"] = pandas.DataFrame(
                    read_features(text_stream),
                    columns=["gene_id", "gene_symbols"],
                )
            elif member_path.suffix == ".mtx":
                result["matrix"] = mmread(text_stream)
            else:
                logger.warning("Unrecognized archive member {}".format(member.name))
    return result


def read_mex_archive_as_anndata(filename=None, fileobj=None):
    """Reads a mex archive into an AnnData structure

    Parameters
    ----------
    filename : Name of file to read
    fileobj : file-like object to read

    one and only one of filename or fileobj needs to be provided.

    Returns
    -------

    AnnData object
    """
    from anndata import AnnData

    result = read_mex_archive(filename=filename, fileobj=fileobj)
    result["metadata"] = remove_manifest_md5s(result["metadata"])
    adata = AnnData(X=result["matrix"].T.tocsc(), uns=result["metadata"])
    adata.obs_names = result["barcodes"]
    adata.var_names = result["features"]["gene_id"]
    adata.var["gene_symbols"] = result["features"]["gene_symbols"].to_numpy()
    return adata


def main(cmdline=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--url", help="url to read from")
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("filename", nargs="?")
    args = parser.parse_args(cmdline)

    if args.url is None and args.filename is None:
        parser.error("Please set a source, either a url or local filename")
    elif args.url is not None and args.filename is not None:
        parser.error("Please only use one source")

    fileobj = None
    if args.url:
        import requests

        req = requests.get(args.url, stream=True)
        req.raise_for_status()
        fileobj = req.raw

    adata = read_mex_archive_as_anndata(filename=args.filename, fileobj=fileobj)
    output = Path(args.output)
    if output.suffix == ".h5ad":
        adata.write_h5ad(args.output)
    elif output.suffix == ".loom":
        adata.write_loom(args.output)
    elif output.suffix == ".zarr":
        adata.write_zarr(args.output)
    else:
        parser.error("Unrecognized output extension")


if __name__ == "__main__":
    main()
