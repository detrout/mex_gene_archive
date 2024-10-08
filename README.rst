|ci-test|

Introduction
============

mex_gene_archive is a minimal file format designed to meet the needs
of archiving sparse gene matrices in a format compatible with the
ENCODE 4 Data Coordination Center.

We had the requirement that a data type result needed to be a single
file and unfortunately the common output format for alignment programs
of the `matrix market exchange`_ use three files. One to store the
coordinates and values of the non-zero sparse matrix elements, one for
the row labels, and one for the column labels.

For IGVF the contents were extended to support any number of matrixes
in one tar file.

The archives include a manifest.tsv file that provides some information
about what is in the archive.

Here is an example subset of a recent manifest file. It contains a list of name value pairs.
In the case of a list metadata value, the additional list items are written as new entries
with a blank value name.

The end of the file will contain a list of the other included files and their md5sum.

One of the weaknesses of tar.gz is the lack of a file content index.

+-----------------------------------------------------+---------------------------------------+
|name                                                 |value                                  |
+-----------------------------------------------------+---------------------------------------+
|type                                                 |MexGeneArchive                         |
|software                                             |kallisto                               |
|software_version                                     |0.50.1                                 |
|input_file_set                                       |TSTDS34582101                          |
|                                                     |TSTDS07432728                          |
|counts_unfiltered_modified/cells_x_genes.barcodes.txt|md5sum:af2e5a8e0f6d2da661e99c64a9656224|
+-----------------------------------------------------+---------------------------------------+

Usage
=====

Reading an archive
------------------

The archive format is fairly simple and started with just archiving
the key matrix market files from a STAR Solo.out directory, with a
simple manifest.tsv file included to help tell different files
apart.

Probably the more useful function is the one that will read an archive
into an anndata structure with the gene features going across the
columns and the cell barcode observations going down across the rows.

.. code-block:: python

   from mex_gene_archive.reader import read_mex_archive_as_anndata

   adata = read_mex_archive_as_anndata("archive.tar.gz")

   req = requests.get(
       "https://www.encodeproject.org/files/ENCFFexample/@@download/ENCFFexample.fastq.gz",
       stream=True)
   adata = read_mex_archive_as_anndata(fileobj=req.raw)


The reader module can also convert archives to anndata directly from
the command line

.. code-block:: bash

   python -m mex_gene_archive.reader -o archive.h5ad archive.tar.gz

   python -m mex_gene_archive.reader -o archive.h5ad \
     --url https://www.encodeproject.org/files/ENCFFexample/@@download/ENCFFexample.fastq.gz


Generating an STAR archive
--------------------------

Possibly you might want to generate an archive file currently only
STAR is directly supported. See archive_star_solo for the full list of
arguments.

.. code-block:: python

   from mex_gene_archive.starsolo import archive_star_solo

   config = {
      "experiment_accession": "ENCSR724KET",
      "description": "snRNA on human adrenal gland.",
      "library_accession": "ENCLB002DZK",
   }
   archive_star_solo("experiment", config)

.. _`matrix market exchange`: https://math.nist.gov/MatrixMarket/
.. |ci-test| image:: https://github.com/detrout/mex_gene_archive/actions/workflows/ci-test.yml/badge.svg
