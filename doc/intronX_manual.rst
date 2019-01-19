==================
Manual for IntronX
==================

IntronX is a program included in pySeqLib_ package to fetch sequence related 
features for introns. 

.. _pySeqLib: https://github.com/huangyh09/pyseqlib


Installation
============

The installation is the same for pySeqLib_. Simply, you download this package 
from this GitHub repository, and unzip it and run the following command line 

::

    python setup.py install

If you don't have the root permission, add ``--user``. When sucessfully 
installed pySeqLib_, you could run ``intronX`` on from your terminal directly. 
Also you can find the codes for this program here_. 
Rquired libraries: ``pysam``,  ``numpy``, ``Cython``.

.. _pySeqLib: https://github.com/huangyh09/pyseqlib
.. _here: https://github.com/huangyh09/pyseqlib/blob/master/pyseqlib/intronX/intronX.py


Input and Ouput
===============

This program is to fetch sequence related features for introns. Therefore, the 
necessary input files are the information of a list of introns, and the genome 
sequence reference file. The genome reference file is in the sandard FASTA file 
format, and the file for intron information should be a plain text file with 
tab separated, i.e., tsv, and should contain 6 or 7 columns, as follows,

  * column 1: gene id
  * column 2: intron id
  * column 3: chromosome id
  * column 4: strand
  * column 5: start site (smaller number)
  * column 6: stop site (larger number)
  * column 7 (optional): branch point site

Then you could input these two files and the directory for outpu, and run it 
like this:

::

  intronX -i $intron_table -f $fasta -o $out_dir

If you don't put an output directory, it will use the one for the genome 
sequence file. By default, you will see a folder ``intronX`` in your output 
directory, and there will be a final result file ``intronX.txt`` and weblogo 
figure ``3ss_seq.pdf``, ``5ss_seq.pdf`` and ``BPs_seq.pdf`` (optional), and 
also a subfolder ``seq`` containing the following files: 

  * for splice site motif: ``3ss_seq.fa``, ``5ss_seq.fa``, ``BPs_seq.fa`` 
    (optional)
  * for k-mers and second structure: ``intron_seq.fa``, ``5ss_BPs_seq.fa`` 
    (optional), ``BPs_3ss_seq.fa`` (optional)
  * for RNAfold second structure results: ``intron_seq.RNAfold.txt``, 
    ``5ss_BPs_seq.RNAfold.txt`` (optional), ``BPs_3ss_seq.RNAfold.txt`` 
    (optional)

As can be seen above, there are 4 types of features 1) length, 2) splice site 
motif strength, 3) second structure energy 4) frequecy of k-mers. In addition 
it could produce the motif logo. 

Here, the second structure is predicted by RNAfold_ and motif logo is produced 
by webLogo_ v3.0. ``RNAfold`` is from the ViennaRNA_ Package, and needed be 
installed beforehand, and available in the `$PATH` variable (see this tutorial_)
Additionally, you could turn off any of these functions off by input 
``--no-RNAfold`` and / or ``--no-weblogo``, respectively. The according output 
files will not be available, consequently.

In addition, you could choose the range of k-mers, by setting 
``--kmer-range 2 4`` (for kmin=2 and kmax=4).

.. _webLogo: https://github.com/WebLogo/weblogo
.. _RNAfold: https://www.tbi.univie.ac.at/RNA/documentation.html
.. _ViennaRNA: https://www.tbi.univie.ac.at/RNA/index.html
.. _tutorial: https://robots.thoughtbot.com/the-magic-behind-configure-make-make-install

There are more parameters for setting (``intronX -h`` always give the version 
you are using)

.. code-block:: html

  Usage: intronX [options]

  Options:
    -h, --help            show this help message and exit
    -i INTRON_TABLE, --intron_table=INTRON_TABLE
                          The intron table file for intron information.
    -f FASTA_FILE, --fasta=FASTA_FILE
                          The fasta file of genome sequence.
    -o OUT_DIR, --out_dir=OUT_DIR
                          The directory for output [default: $fasta].

    Optional arguments:
      --kmer-range=KMER_RANGE
                          The min and max K in k-mers. [default: 1 3]
      --no-RNAfold        No second strucutre for intron sequences
      --no-weblogo        No weblogo figures for motifs


Examples
========

There are some examples available here: 
https://sourceforge.net/projects/pyseqlib/files/intronX-example/

- Example to for introns in yeast with or without branch point (bash code and 
  data): intronYeast.zip_

.. _intronYeast.zip: http://ufpr.dl.sourceforge.net/project/pyseqlib/intronX-example/intronYeast.zip