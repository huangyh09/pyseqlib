# This file is designed to extract intron sequence, 5 splice site,
# 3 splice site and branch point (if given at the 7th column).
# The input is the intron_table, which contains 7 columns or 6 (if
# no branch points): Name, chromasome, strand, start, stop, bp.


import os
import sys
import numpy as np
from optparse import OptionParser
from pyseqlib import FastaFile, rev_seq, cDNA2RNA, fasta_write

if __name__ == "__main__":
    print("Welcome to IntronSeq!")

    # 0. parse command line options
    parser = OptionParser()
    parser.add_option("--intron_table", "-i", dest="intron_table", default=None,
        help="The intron table file for intron information.")
    parser.add_option("--fasta", "-f", dest="fasta_file", default=None,
        help="The fasta file of genome sequence.")
    parser.add_option("--out_dir", "-o", dest="out_dir",  
        default=None, help="The directory for output [default: $gff_dir].")

    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("use -h or --help for help on argument.")
        sys.exit(1)

    if options.intron_table is None:
        print("Error: need --intron_table for intron information.")
        sys.exit(1)
    else:
        intron_info = np.loadtxt(options.intron_table, delimiter="\t", 
            dtype="str", skiprows=1)

    if intron_info.shape[1] == 5: bp_exist = False
    else: bp_exist = True

    if options.fasta_file is None:
        print("Error: need --fasta for fasta file of genome sequence.")
        sys.exit(1)
    else:
        fastaFile = FastaFile(options.fasta_file)

    if options.out_dir is None:
        out_dir = os.path.dirname(options.fasta_file) + "/intronSeq"
    else:
        out_dir = options.out_dir #+ "/intronSeq"
    try:
        os.stat(out_dir)
    except:
        os.mkdir(out_dir)


    fid1 = open(out_dir + "/intron.fa", "w")
    fid2 = open(out_dir + "/5ss.fa", "w")
    fid3 = open(out_dir + "/3ss.fa", "w")
    if bp_exist: fid4 = open(out_dir + "/bp.fa", "w")

    for i in range(intron_info.shape[0]):
        chrom = str(intron_info[i,1])
        strand = str(intron_info[i,2])
        UPs = int(intron_info[i,3])
        DNs = int(intron_info[i,4])
        if bp_exist: BPs = int(intron_info[i,5])

        if strand == "+" or strand == "1":
            _seq_intron = fastaFile.get_seq(chrom, UPs, DNs)
            _seq_5ss = fastaFile.get_seq(chrom, UPs-4,  UPs+7)
            _seq_3ss = fastaFile.get_seq(chrom, DNs-16, DNs+4)
            if bp_exist: _seq_BPs = fastaFile.get_seq(chrom, BPs-7,  BPs+3)
        else :
            _seq_intron = rev_seq(fastaFile.get_seq(chrom, UPs, DNs))
            _seq_5ss = rev_seq(fastaFile.get_seq(chrom, DNs-7,  DNs+4))
            _seq_3ss = rev_seq(fastaFile.get_seq(chrom, UPs-4, UPs+16))
            if bp_exist: _seq_BPs = rev_seq(fastaFile.get_seq(chrom, BPs-3,  BPs+7))
        
        # fasta_write(fid1, cDNA2RNA(_seq_intron),  intron_info[i,0], length=60)
        fasta_write(fid1, _seq_intron,  intron_info[i,0], length=60)
        fasta_write(fid2, _seq_5ss, intron_info[i,0], length=60)
        fasta_write(fid3, _seq_3ss, intron_info[i,0], length=60)
        if bp_exist: fasta_write(fid4, _seq_BPs, intron_info[i,0], length=60)

    fid1.close()
    fid2.close()
    fid3.close()
    if bp_exist: fid4.close()

