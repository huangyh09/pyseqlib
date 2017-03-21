# This file is to randomly sample reads from fastq file.

import io
import os
import sys
import gzip
import time
import numpy as np
from optparse import OptionParser

def main():
    # parse command line options
    parser = OptionParser()
    parser.add_option("--inFile", "-i", dest="input_file", default=None,
        help="Full path of input fastq file.")
    parser.add_option("--outFile", "-o", dest="out_file", default=None, 
        help="Full path of output file.")
    parser.add_option("--numSample", "-n", dest="num_sample", default=None, 
        help="Number of sampled reads.")
    parser.add_option("--ungzip", dest="num_sample", default=None, 
        help="Number of sampled reads.")
    parser.add_option("--plainOut", action="store_true", dest="plain_out", 
        default=False, help="Save plain text file for output, not gzip.")

    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("Welcome to Fastq-Sample!\n")
        print("use -h or --help for help on argument.")
        sys.exit(1)

    if options.input_file is None:
        print("Error: need input fastq file.")
        sys.exit(1)
    else:
        input_file = options.input_file
        file_name  = ".".join(os.path.basename(input_file).split(".")[:-1])

    if options.num_sample is None:
        print("Error: need -n numSample for number of sampled reads.")
        sys.exit(1)
    else:
        num_sample = int(options.num_sample)

    if options.out_file is None:
        out_file = os.path.dirname(input_file) + "/random%d.fq.gz" %(num_sample)
    else:
        out_file = options.out_file

    START_TIME = time.time()

    # counting total reads
    total_reads = 0
    try:
        infile = gzip.open(input_file, 'rb')
        with io.BufferedReader(infile) as f:
            for line in f:
                total_reads += 1
        ftype = "gzip"
    except:
        infile = open(input_file, 'rb')
        for line in infile:
            total_reads += 1
        ftype = "plain"
    infile.close()
    total_reads = total_reads / 4
    sys.stdout.write('\r[Fastq-Sample] Sample %d out of %d reads.' 
        %(num_sample, total_reads))
    sys.stdout.flush()

    # generate random reads index
    idx_out = np.random.permutation(total_reads)[:num_sample]
    idx_out = np.sort(idx_out)
    
    ## print idx_out, len(idx_out)
    ## idx_out = np.arange(100)

    # output sampled reads
    if ftype == "gzip":
        infile = io.BufferedReader(gzip.open(input_file, 'rb'))
    else:
        infile = open(input_file, 'rb')

    if options.plain_out is True: 
        outfile = open(out_file, "w")
    else:
        outfile = gzip.open(out_file, "w")

    outCNT, lineCNT = 0, -1
    for line in infile:
        lineCNT += 1
        if int(lineCNT / 4) == idx_out[outCNT]:
            outfile.writelines(line)
            if lineCNT % 4 == 3: outCNT += 1
        if outCNT >= len(idx_out):
            break
    infile.close()
    outfile.close()

    run_time = time.time() - START_TIME
    sys.stdout.write(('\r[Fastq-Sample] Sample %d out of %d reads. '
        'Done in %.2f sec.' %(num_sample, total_reads, run_time)))
    sys.stdout.flush()
    print("")


if __name__ == "__main__":
    main()