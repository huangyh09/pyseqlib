pySeqLib
========

A small python library to process sequencing data. It is under active development.


What's inside?
---------------

* ``intronX``: a program to fetch sequence related features for introns: 
  1) length, 2) splice site motif strength, 3) second structure energy 
  4) frequecy of k-mers. See manual_ and example_.

* ``pymfold``: a python wrap of mfold for calculate energy of RNA secondary 
  structure. We got installation issue from mfold, similar as this_. Therefore,
  we recommend using `pyRNAfold` instead.

* ``pyRNAfold``: a python wrap of for another predictor, RNAfold, of RNA 
  secondary structure energy. It is from the ViennaRNA_ package.

* ``motif-score``: a motif score calculate

* a sequence mapper to find lariat in RNA-seq reads (under test)

* a sam file resampling method

.. _this: http://unafold.rna.albany.edu/?q=node/927
.. _ViennaRNA: https://www.tbi.univie.ac.at/RNA/index.html
.. _manual: https://github.com/huangyh09/pyseqlib/blob/master/doc/intronX_manual.rst
.. _example: https://sourceforge.net/projects/pyseqlib/files/intronX-example/


How to install?
---------------

Pyseqlib was initially developed in Python 2 environment, hence best to be used 
in Py2 environment. By using Anaconda platform, no matter Py2 or Py 3, it is 
easy to set up a conda environment with Py2, for example by following commond:

.. code-block:: bash

   conda create -n Py2 python=2.7 numpy==1.15.4 pysam==0.15.2 Cython==0.29.3 matplotlib==2.2.3
   source activate Py2

Once you are in a Python 2 environment, you can download the codes from this 
github repository and then run the following command line:

::

    python setup.py install

If you don't have the root permission, add ``--user``.

Rquired libraries: ``pysam``,  ``numpy``, ``Cython``.

