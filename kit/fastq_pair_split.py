# This file is to split the paired-end reads in a single fastq 
# file generated from FluxSimulator.

import os
import sys
from optparse import OptionParser

def main():
    print("Welcome to Fastq-Pair-Split!")

    # parse command line options
    parser = OptionParser()
    parser.add_option("--input", "-i", dest="input", default=None,
        help="The input fastq file, containing both mates.")
    parser.add_option("--outDir", "-o", dest="out_dir", default=None, 
        help="The directory for output [default: $input].")

    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("use -h or --help for help on argument.")
        sys.exit(1)

    if options.input is None:
        print("Error: need input fastq file.")
        sys.exit(1)
    else:
        fastq = options.input
        fname = ".".join(os.path.basename(fastq).split(".")[:-1])

    if options.out_dir is None:
        out_dir = os.path.dirname(fastq) + "/"
    else:
        out_dir = options.out_dir + "/"
    try:
        os.stat(out_dir)
    except:
        os.mkdir(out_dir)

    fid1 = open(out_dir + fname + "_1.fq", "w")
    fid2 = open(out_dir + fname + "_2.fq", "w")

    # process the data
    cnt = 0
    read = []
    read_num = 1
    with open(fastq, "r") as infile:
        for line in infile:
            cnt += 1
            read.append(line.rstrip())
            if cnt % 4 == 0:
                if read[0][-1] == "1":
                    fid1.writelines("@FluxSimulator.read:%d_1\n" %read_num)
                    fid1.writelines(read[1] + "\n")
                    fid1.writelines("+%s\n" %read[0][1:])
                    fid1.writelines(read[3] + "\n")
                else:
                    fid2.writelines("@FluxSimulator.read:%d_2\n" %read_num)
                    fid2.writelines(read[1] + "\n")
                    fid2.writelines("+%s\n" %read[0][1:])
                    fid2.writelines(read[3] + "\n")
                read = []
            if cnt % 8 == 0:
                read_num += 1
    fid1.close()
    fid2.close()


if __name__ == "__main__":
    main()