Introduction
============

`mex_gene_archive` is a small library to create tar archives of matrix market
files with some minimal metadata embedded in the archive.

The project was written because the ENCODE DCC had a requirement that
a single result be contained in a single file, but unfortunately the
sparse matrix format `Matrix Market Exchange`_ uses three text files
to define a matrix.

A column label text file, a row label text file, and file encoding
sparse matrix as a set of coordinates followed by a value.

The full matrix market specification documents variant options that
can further reduce the number of values present for matrices that have
certain types of symmetries.

_`Matrix Market`: https://math.nist.gov/MatrixMarket/
