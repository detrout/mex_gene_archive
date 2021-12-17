Changelog
=========

Release 0.2.0
-------------

This version is much better at removing all the possible sources of
variation in the tar and gzip file formats so if the input files are
the same the file hashes should also match.

I used the older way of initialzing the numpy random number
generator to support numpy as far back as as 1.16.

Fixed a problem wehre scipy's mmread expected a bytestream and not a
textstream
