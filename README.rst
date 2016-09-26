pySeqLib
========

A small python library to process sequencing data. It is under active development.


What's inside?
---------------

* a python wrap of mfold for calculate energy of RNA secondary structure.

* a motif score calculate

* a sequence mapper to find lariat in RNA-seq reads (under test)

* a sam file resampling method


How to install?
---------------

Download the codes from this github repository and then run the following command line:

::

    python setup.py install

If you don't have the root permission, add ``--user``.

Rquired libraries: ``pysam``,  ``numpy``, ``Cython``.

