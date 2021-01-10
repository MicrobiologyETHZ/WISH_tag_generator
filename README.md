# tnseq_tags

A program to generate tags for both QPCR and HTS.

The method is approximately:

* Generate random sequence
* Check GC content 45% - 55%
* Check for palindrome
* Align to database, allowing 2 errors
* Check for mispriming
* Score with Primer3

The method does not:

* Check the tags against one another
* Check the tags against Illumina barcodes/primers
