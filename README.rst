======
idemux
======
.. image:: https://badge.fury.io/py/idemux.svg
   :target: https://badge.fury.io/py/idemux
   :alt: Latest Version

.. image:: https://travis-ci.org/lexogen-tools/idemux.svg?branch=master
   :target: https://travis-ci.org/lexogen-tools/idemux

.. image:: https://coveralls.io/repos/github/lexogen-tools/idemux/badge.svg?branch=master
   :target: https://coveralls.io/github/lexogen-tools/idemux?branch=master

Idemux is a command line tool designed to demultiplex paired-end fastq files from 
`QuantSeq-Pool <https://www.lexogen.com/quantseq-pool-sample-barcoded-3mrna-sequencing/>`_.

Idemux can demultiplex based on i7, i5 and i1 inline barcodes. It has a built in 
sequencing error correction for `Lexogen indices <https://www.lexogen.com/indexing/12nt-dual-indexing-kits/>`_, 
allowing you to retain more reads from your sequencing experiment


Installation
-----------

Idemux is available on pypi. To install idemux simply install with pip.

``pip install idemux``

Idemux will soon also be available via bioconda 


How-to
-------


In order to run idemux you 




Features
--------

* TODO

Help
------

* OSError: [Errno 24] Too many open files
* ulimit -n number_of_samples_x2 (in 1024 increments)

