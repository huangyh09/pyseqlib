# This program is to return motif scores, and save motif logo

# Author: Yuanhua Huang, Y.Huang@ed.ac.uk
# Date: 2015-09-19

import os
import sys
import numpy as np
from optparse import OptionParser

import pyximport; pyximport.install()
from .utils.fasta import LightFasta, motif_score

def main():
    print "Welcome to motif_score."

    #part 0. parse command line options
    parser = OptionParser()
    parser.add_option("--motif_file", "-f", dest="motif_file",
        help="The motif file in fasta or msa.")
    parser.add_option("--pwm_file", "-w", dest="pwm_file", default=None,
        help="The input position weights matrix, in fasta file.")
    parser.add_option("--motif_file_type", dest="motif_file_type", default="fasta",
        help="The type of motif file, fasta or msa: [default: %default].")
    parser.add_option("--out_file", "-o", dest="out_file",
        default="motif_score_out.txt", help="The out file results.")

    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("use -h or --help for help on argument.")
        sys.exit(1)

    fasta = LightFasta(options.motif_file)
    if options.pwm_file is not None:
        pwm_seq = LightFasta(options.pwm_file).seq
    else:
        pwm_seq = None
    out_file = options.out_file

    scores = motif_score(fasta.seq, pwm_seq)
    fid = open(out_file, "w")
    fid.writelines("id\tscore\n")
    for i in range(len(fasta.ref)):
        fid.writelines("%s\t%.2f\n" %(fasta.ref[i], scores[i]))
    fid.close()


if __name__ == "__main__":
    main()
