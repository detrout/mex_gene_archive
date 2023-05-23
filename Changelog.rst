
Changelog
=========

Release 0.2.3
-------------

The code to generate the md5sums for the merged matrix file didn't
call .hexdigest() so the md5sum for the merged sparse matrix file was
the string rerepresentation of the python md5 object.


Release 0.2.2
-------------

Add function to write an AnnData archive to a mex_gene_archive file.

def write_anndata_as_mex_archive(
    adata,
    quantification="GeneFull",
    multiread="Unique",
    matrix="raw",
    *,
    destination=None)

Release 0.2.1
-------------

The manifest.tsv file is now being written into the same logical
directory as the matrix being archived, to make it easier to extract
multiple matrices.

read_mex_archive_as_anndata was adjusted to treat the splice junction
matrix (SJ_Unique_raw..tar.gz) feature table differently from the gene
features. Now it returns all the information about the splice
junctions that STAR provides.

Internally there was some reorganization to support archiving matrices
from memory instead of disk for a feature that turned out to be better
implemented elsewhere.

Release 0.2.0
-------------

This version is much better at removing all the possible sources of
variation in the tar and gzip file formats so if the input files are
the same the file hashes should also match.

I used the older way of initialzing the numpy random number
generator to support numpy as far back as as 1.16.

Fixed a problem wehre scipy's mmread expected a bytestream and not a
textstream
