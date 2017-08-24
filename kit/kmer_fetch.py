# This program is to fetch k-mers from a sequence file in fasta formate.
# Note: Please make sure the sequence is in the same order, especially
# minus strand, been reversed properly.

# Author: Yuanhua Huang, Y.Huang@ed.ac.uk
# Date: 2017-08-23

import os
import sys
import time
import numpy as np
import multiprocessing
from optparse import OptionParser, OptionGroup

import pyximport; pyximport.install()
from .utils.fasta import LightFasta, get_kmer_all, get_motif, rev_seq


PROCESSED = 0
TOTAL_GENE = 0
START_TIME = time.time()
    
def show_progress(RV=None):
    global PROCESSED, TOTAL_GENE, START_TIME
    if RV is not None: 
        PROCESSED += 1
        bar_len = 20
        run_time = time.time() - START_TIME
        percents = 100.0 * PROCESSED / TOTAL_GENE
        filled_len = int(bar_len * percents / 100)
        bar = '=' * filled_len + '-' * (bar_len - filled_len)
        
        sys.stdout.write('\r[kmer-fetch] [%s] %.1f%% done in %.1f sec.' 
            % (bar, percents, run_time))
        sys.stdout.flush()
    return RV


def main():
    print "Welcome to kmer-fetch."

    #part 0. parse command line options
    parser = OptionParser()
    parser.add_option("--fasta", "-f", dest="fasta_file", default=None,
        help="The sequence file in fasta format.")
    parser.add_option("--out_file", "-o", dest="out_file", default=None, 
        help="The output file for kmer counts in tsv format.")

    group = OptionGroup(parser, "Optional arguments")
    group.add_option("--nproc", "-p", type="int", dest="nproc", default="1",
        help="Number of subprocesses [default: %default]")
    group.add_option("--kmin", type="int", dest="kmin", default="1",
        help="Minimum of k in k-mers [default: %default]")
    group.add_option("--kmax", type="int", dest="kmax", default="3",
        help="Maximum of k in k-mers [default: %default]")
    group.add_option("--kmer_seq", dest="kmer_seq", default="ATGC",
        help=("The list of letters in the sequences."))
    group.add_option("--kmer_mode", dest="kmer_mode", default="frequency",
        help=("The mode for the kmer counting: frequency, count, normalized."))
    parser.add_option_group(group)

    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("use -h or --help for help on argument.")
        sys.exit(1)

    if options.fasta_file is None:
        print("Error: need --fasta for fasta file of candidate sequences.")
        sys.exit(1)
    else:
        fastaFile = LightFasta(options.fasta_file)
        _ref, _seq = fastaFile.ref, fastaFile.seq

    if options.out_file is None:
        out_file = os.path.dirname(options.fasta_file) + "/kmer_counts.tsv"
    else:
        out_file = options.out_file #+ "/intronSeq"
    try:
        os.stat(os.path.dirname(out_file))
    except:
        os.mkdir(os.path.dirname(out_file))

    kmin = options.kmin
    kmax = options.kmax
    kmer_seq = options.kmer_seq
    kmer_lst = get_kmer_all(kmax=kmax, kmin=kmin, seqs=kmer_seq)
    kmer_frq = np.zeros((len(_ref), len(kmer_lst)))
    kmer_mode = options.kmer_mode

    # count k-mer
    global TOTAL_GENE
    TOTAL_GENE = len(kmer_lst) * len(_ref)
    if options.nproc <= 1:
        for i in range(len(_ref)):
            for j in range(len(kmer_lst)):
                kmer_frq[i,j] = get_motif(_seq[i], kmer_lst[j], mode=kmer_mode)
                show_progress(kmer_frq[i,j])
    else:
        pool = multiprocessing.Pool(processes=options.nproc)
        result = []
        for i in range(len(_ref)):
            for j in range(len(kmer_lst)):
                result.append(pool.apply_async(get_motif, (_seq[i], kmer_lst[j], 
                    kmer_mode), callback=show_progress))
        pool.close()
        pool.join()
        cnt = -1
        for res in result:
            cnt += 1
            n, m = int(cnt/len(kmer_lst)), cnt%len(kmer_lst)
            kmer_frq[n,m] = res.get()

    # save output
    fid = open(out_file, "w")
    headline = "\t".join(["seq_id"] + kmer_lst)
    fid.writelines(headline + "\n")
    for i in range(len(_ref)):
        aline = _ref[i] + "\t"
        aline += "\t".join(["%.2e" %x for x in kmer_frq[i,:]])
        fid.writelines(aline + "\n")
    fid.close()
    print("")


if __name__ == "__main__":
    main()
