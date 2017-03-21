# dump genome.fa sequence into chrom_xx.fa

import os
import sys
from optparse import OptionParser

def main():
    print("Welcome to Genome-dump!")

    # parse command line options
    parser = OptionParser()
    parser.add_option("--genome", "-g", dest="genome", default=None,
        help="The genome file in fasta.")
    parser.add_option("--outDir", "-o", dest="out_dir", default=None, 
        help="The directory for output [default: $genome_dir/chroms].")
    parser.add_option("--lineLen", "-l", dest="lineLen", default="50", 
        help="The length of each line [default: %default].")

    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("use -h or --help for help on argument.")
        sys.exit(1)

    if options.genome is None:
        print("Error: need genome sequence in fasta.")
        sys.exit(1)
    else:
        genome = options.genome

    if options.out_dir is None:
        out_dir = os.path.dirname(genome) + "/chroms"
    else:
        out_dir = options.out_dir
    try:
        os.stat(out_dir)
    except:
        os.mkdir(out_dir)

    lineLen = int(options.lineLen)

    # process the data
    chrom_begin = False
    with open(genome, "r") as infile:
        for line in infile:
            if line.startswith(">"):
                prevLen = 0
                if chrom_begin is True:
                    fid.close()
                else:
                    chrom_begin = True
                fid = open(out_dir + line.split(" ")[0][1:] + ".fa", "w")
                fid.writelines(line.split(" ")[0] + "\n")

            elif chrom_begin is True:
                line = line.rstrip()
                chromLen = len(line)
                curr_len = lineLen - prevLen
                if curr_len >= chromLen:
                    curr_loc = chromLen
                    fid.writelines(line)
                    prevLen = prevLen + curr_len
                else:
                    fid.writelines(line[:curr_len] + "\n")
                    curr_loc = curr_len
                    while curr_loc < chromLen-lineLen:
                        fid.writelines(line[curr_loc : curr_loc+lineLen] + "\n")
                        curr_loc = curr_loc + lineLen
                    prevLen = chromLen - curr_loc
                    fid.writelines(line[curr_loc : curr_loc+prevLen])
    fid.close()


if __name__ == "__main__":
    main()