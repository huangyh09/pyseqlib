# This file is developed for convert fasta file(s) into new fasta file(s)
# that have defined length of each line. This program has tow options: 
# (1) combine many fasta file into a single fasta file;
# (2) split a big fasta file into many small fasta files with a single
# ref_id in each file. For this case, --split is required.

# Example: python fasta_split.py -i chr1.fa,chr2.fa -o genome.fa


import os
import sys
from optparse import OptionParser

def main():
    # parse command line options
    parser = OptionParser()
    parser.add_option("--inFasta", "-i", dest="in_files", default=None,
        help="Input fasta file(s), e.g., C1.fa,C2.fa")
    parser.add_option("--outDir", "-o", dest="out_dir", default=None, 
        help="Full path of output directory [default: $inFasta[0]/faConv].")
    parser.add_option("--lineLen", "-l", dest="line_len", default="50", 
        help="Length of each line [default: %default].")
    parser.add_option("--split", action="store_true", dest="is_split",
        default=False, help="Split the output into many fasta files "
        "with each ref_id")

    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("Welcome to Fasta-Convert!\n")
        print("use -h or --help for help on argument.")
        sys.exit(1)

    if options.in_files is None:
        print("Error: need input fasta file.")
        sys.exit(1)
    else:
        in_files = options.in_files.split(",")
        print("Converting %d fasta files... " %(len(in_files)))

    is_split = options.is_split
    line_len = int(options.line_len)

    if options.out_dir is None:
        if is_split is True:
            out_dir = os.path.dirname(in_files[0]) + "/faConv"
        else:
            out_dir = os.path.dirname(os.path.abspath(in_files[0]))
    else:
        out_dir = options.out_dir
    try:
        os.stat(out_dir)
    except:
        os.mkdir(out_dir)

    # process the data
    chrom_begin = False
    for aFasta in in_files:
        with open(aFasta, "r") as infile:
            for line in infile:
                if line.startswith(">"):
                    prevLen = 0
                    if is_split is False:
                        if chrom_begin is False:
                            fid = open(out_dir + "/converted.fa", "w")
                            chrom_begin = True
                        else:
                            fid.writelines("\n")
                        fid.writelines(line)
                    else:
                        if chrom_begin is True:
                            fid.close()
                        else:
                            chrom_begin = True
                        fid = open(out_dir + line.split(" ")[0][1:] + ".fa", "w")
                        fid.writelines(line.split(" ")[0] + "\n")

                elif chrom_begin is True:
                    line = line.rstrip()
                    chromLen = len(line)
                    curr_len = line_len - prevLen
                    if curr_len >= chromLen:
                        curr_loc = chromLen
                        fid.writelines(line)
                        prevLen = prevLen + curr_len
                    else:
                        fid.writelines(line[:curr_len] + "\n")
                        curr_loc = curr_len
                        while curr_loc < chromLen-line_len:
                            fid.writelines(line[curr_loc : curr_loc+line_len] + "\n")
                            curr_loc = curr_loc + line_len
                        prevLen = chromLen - curr_loc
                        fid.writelines(line[curr_loc : curr_loc+prevLen])
    fid.close()
    print("Converted!")


if __name__ == "__main__":
    main()