pySeqLib
========

A small python library to process sequencing data. It is under active development.


What's inside?
---------------

* ``intronX``: a program to fetch sequence related features for introns: 
  1) length, 2) splice site motif strength, 3) second structure energy 
  4) frequecy of k-mers. See manual_ and example_.

* ``pymfold``: a python wrap of mfold for calculate energy of RNA secondary 
  structure.

* ``motif-score``: a motif score calculate

* a sequence mapper to find lariat in RNA-seq reads (under test)

* a sam file resampling method


.. _manual: https://github.com/huangyh09/pyseqlib/blob/master/doc/intronX_manual.rst
.. _example: https://sourceforge.net/projects/pyseqlib/files/intronX-example/


How to install?
---------------

Download the codes from this github repository and then run the following command line:

::

    python setup.py install

If you don't have the root permission, add ``--user``.

Rquired libraries: ``pysam``,  ``numpy``, ``Cython``.

