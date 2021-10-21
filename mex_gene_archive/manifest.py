import csv
import hashlib

#######
# Functions for making manifest
#
IO_BUFFER = 2 ** 10


def compute_md5sums(filenames):
    results = []
    for f in filenames:
        with open(f, "rb") as instream:
            md5 = compute_stream_md5sum(instream)
        results.append((f, md5.hexdigest()))
    return results


def compute_stream_md5sum(instream):
    """Compute the md5sum of a file-like object
    """
    md5 = hashlib.md5()
    readable = instream.readable()
    while readable:
        read_block = instream.read(IO_BUFFER)
        if len(read_block) == 0:
            readable = False
        else:
            md5.update(read_block)
    return md5


def create_metadata(config, md5s):
    metadata = {
        "type": "MexGeneArchive_v1",
        "output_type": config['output_type'],
        "experiment_accession": config.get("experiment_accession", ""),
        "description": config.get("description", ""),
        "library_accession": config.get("library_accession", ""),
        "analysis_version": config.get("analysis_version", ""),
    }

    for filename, md5 in md5s:
        metadata[filename] = "md5sum:{}".format(md5)

    return metadata


def write_manifest(outstream, config):
    writer = csv.writer(outstream, delimiter="\t")
    writer.writerow(["name", "value"])
    for key in config:
        writer.writerow([key, config[key]])
    return outstream


def read_manifest(instream):
    reader = csv.reader(instream, delimiter="\t")
    header = None
    metadata = {}
    for row in reader:
        if header is None:
            header = row
        else:
            metadata[row[0]] = row[1]
    return metadata
