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

.. _`matrix market exchange`: https://math.nist.gov/MatrixMarket/
