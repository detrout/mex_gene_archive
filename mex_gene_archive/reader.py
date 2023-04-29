import argparse
from io import TextIOWrapper
import logging
import numpy
from pathlib import Path
import pandas
from scipy.io import mmread
import tarfile

from .manifest import read_manifest


logger = logging.getLogger(__name__)


def read_barcodes(instream):
    for line in instream:
        yield line.rstrip()


def read_gene_features(instream):
    info = pandas.read_csv(
        instream,
        sep="\t",
        usecols=[0, 1],
        dtype=str,
        header=None,
        names=["gene_id", "gene_symbols"],
    )

    for i, row in info.iterrows():
        if pandas.isnull(row["gene_symbols"]):
            info["gene_symbols"].iloc[i] = row["gene_id"]
    return info


def read_sj_features(instream):
    return pandas.read_csv(
        instream,
        sep="\t",
        header=None,
        names=[
            "contig",
            "start",
            "end",
            "strand",
            "intron_motif",
            "annotated",
            "unique",
            "multi",
            "overhang",
        ],
    )


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
                if member_path.parts[0] == "SJ":
                    result["features"] = read_sj_features(text_stream)
                else:
                    result["features"] = read_gene_features(text_stream)
            elif member_path.suffix == ".mtx":
                result["matrix"] = mmread(archive.extractfile(member))
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
    if isinstance(result["matrix"], numpy.ndarray):
        X = result["matrix"].T
    else:
        X = result["matrix"].T.tocsc()
    adata = AnnData(X=X, uns=result["metadata"])
    adata.obs_names = result["barcodes"]

    for column in result["features"]:
        if column == "gene_id":
            adata.var_names = result["features"][column].to_numpy()
        else:
            adata.var[column] = result["features"][column].to_numpy()
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
