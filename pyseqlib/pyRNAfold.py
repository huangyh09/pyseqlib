# This program is to run RNAfold with a fasta file and return the minimum free 
# energy values that are predicted by RNAfold, a function from ViennaRNA.
# Please install the ViennaRNA package before using this wrap function!
# https://www.tbi.univie.ac.at/RNA/documentation.html

# Note, when install ViennaRNA, make sure the installed directory is in the 
# $PATH variable. You can set by ./configure --prefix you_PATH

# demo: python pyRNAfold.py -f in_file.fasta -o out_file.txt

# Author: Yuanhua Huang, yuanhua@ebi.ac.uk
# Date: 2019-01-19

import os
import sys
import shutil
import subprocess
from optparse import OptionParser

import pyximport; pyximport.install()
from .utils.fasta import LightFasta

# import os
# DEVNULL = open(os.devnull, 'wb')


def parse_RNAfold_out(RNAfold_out):
    """Parse the RNAfold output stream.
    The format is 
    >RNA_id
    [seqence]
    [structure ([potential_blank]score)]
    """
    if type(RNAfold_out) == bytes:
        RNAfold_out = RNAfold_out.decode("utf-8")
    lines = RNAfold_out.split("\n")
    
    RNA_ids, seq_len, energys  = [], [], []
    for _line in lines:
        if (len(_line) == 0):
            continue
        if _line.startswith(">"):
            RNA_ids.append(_line[1:])
        else:
            _line_list = _line.split(" (")
            if len(_line_list) == 1:
                seq_len.append(len(_line))
            else:
                energys.append(float(_line_list[1][:-1]))
            if (len(RNA_ids) < len(energys)):
                RNA_ids.append("RNA_id%d" %len(energys))

    if len(RNA_ids) != len(energys) or len(RNA_ids) != len(seq_len):
        print("Warning: not matched between of RNA ids, sequences, and energy scores.")
    return RNA_ids, seq_len, energys


def get_RNAfold(fasta_file, out_file=None):
    """get energy score from mfold. Make sure you have installed mfold.
    """
    bashCommand = "RNAfold --noPS %s" %(fasta_file) 
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)        
    RNAfold_output = process.communicate()[0]
    RNA_ids, seq_len, energys =  parse_RNAfold_out(RNAfold_output)
    
    if out_file is None:
        out_file = ".".join(fasta_file.split(".")[:-1]) + ".RNAfold.txt"
    fout = open(out_file, "w")
    fout.writelines("IDs\tlength\tenergy\n")
    for i in range(len(energys)):
        rt_line = [RNA_ids[i], str(seq_len[i]), str(energys[i])]
        fout.writelines("\t".join(rt_line) + "\n")
    fout.close()

    return [str(x) for x in energys]


def main():
    print("Welcome to pyRNAfold, a python wrap for RNAfold from ViennaRNA.")

    #part 0. parse command line options
    parser = OptionParser()
    parser.add_option("--fasta_file", "-f", dest="fasta_file",
        help="The input fasta file for the minimum free energy by mfold.")
    parser.add_option("--out_file", "-o", dest="out_file",
        default="out_file.txt", help="The out file results.")

    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("use -h or --help for help on argument.")
        sys.exit(1)

    out_file = options.out_file
    fasta_file = options.fasta_file

    get_RNAfold(fasta_file, out_file)


if __name__ == "__main__":
    main()
