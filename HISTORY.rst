=======
History
=======


0.1.6 (2021-07-09)
------------------

* Bug fix: Idemux threw an error when both **i7** and **i5** indices were present in the fastq header, but only **i1** demultiplexing should be performed.
* Bug fix: Idemux was not demultiplexing files correctly when both **i7** and **i5** barcodes were present in the fastq header, but only **i5** and **i1** demultiplexing should be performed.
* Enhancement: Added more tests to catch bugs like listed above.


0.1.5 (2020-11-24)
------------------

* Bug fix: Idemux now prints version properly
* Bug fix: README.rst contained some formatting errors
* Bug fix: Broken licence link on pypi now works


0.1.4 (2020-11-12)
------------------

* Bug fix: Demultiplexing with i1 barcodes only raised an incorrect exception (when no barcodes were present in the fastq header)


0.1.3 (2020-08-21)
------------------

* First release on PyPI


0.1.2 (2020-08-21)
------------------

* Bumped version number to avoid upload conflicts


0.1.1 (2020-08-21)
------------------

* Fixed rst files with linter
* First release on test PyPI.


0.1.0 (2020-08-21)
------------------

* First building version.
