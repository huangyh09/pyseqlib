# #
# should output lariat reads!!!


import os
import sys
import time
import numpy as np
import multiprocessing
from optparse import OptionParser

import pyximport; pyximport.install()
from .utils.fasta import RNA2cDNA, rev_seq
from .utils.seq_map import map_lariat_reads


PROCESSED = 0
MAPPED_READ = 0
START_TIME = time.time()

def show_progress(RV=None):
    global PROCESSED, MAPPED_READ, START_TIME, FID
    if len(RV["bp"]) > 0:
        MAPPED_READ += 1
        for bp in RV["bp"]:
            print bp
            FID1.writelines("\t".join(str(x) for x in bp) + "\n")
        for read in RV["read"]:
            FID2.writelines(read)
        if len(RV["read"]) > 0: FID2.writelines("\n")
    
    PROCESSED += 1
    if PROCESSED % 100 == 0:
        run_time = time.time() - START_TIME
        sys.stdout.write("\r%d lariat reads found from %d reads in %.2f s." 
            %(MAPPED_READ, PROCESSED, run_time))
        sys.stdout.flush()
    return RV


def load_fasta(fasta):
    _seq, ref, seq = "", [], []
    with open(fasta, "r") as infile:
        for line in infile:
            line = line.rstrip()
            if len(ref) == 0 and line[0] != ">":
                # comment line
                continue
            if line.startswith(">"):
                ref.append(line[1:])
                if _seq != "": 
                    _seq = RNA2cDNA(_seq)
                    seq.append([_seq, _seq[::-1], 
                             rev_seq(_seq)[::-1], rev_seq(_seq)])
                    _seq = ""
            else:
                _seq = _seq + line
        _seq = RNA2cDNA(_seq)
        seq.append([_seq, _seq[::-1], rev_seq(_seq)[::-1], rev_seq(_seq)])
    return ref, seq


def main():
    print("Welcome to Lariat-Map!")

    # 0. parse command line options
    parser = OptionParser()
    parser.add_option("--fastq", "-r", dest="fastq", default=None,
        help="The fastq file for single-end read.")
    parser.add_option("--intronRef", "-x", dest="intronRef", default=None,
        help="The fasta file for intron sequences.")
    parser.add_option("--outDir", "-o", dest="out_dir", default=None, 
        help="The directory for output [default: same as ref].")

    parser.add_option("--nproc", "-p", dest="nproc", default=1,
        help="The number of processors to use [default: %default].")
    parser.add_option("--overHang", "-H", dest="overhang", default=14,
        help="The minimum overhang [default: %default].")
    parser.add_option("--misMatch", "-e", dest="mismatch", default=3,
        help="The maximum mis-match within overhang [default: %default].")
    parser.add_option("--endCut", "-c", dest="end_cut", default=1,
        help=("The locus to cut near 5ss and bp, in case of different start "
        "index [default: %default]."))
    parser.add_option("--maxMap", "-m", dest="max_map", default=1,
        help="The maximum mapped positions for each read [default: %default].")

    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("use -h or --help for help on argument.")
        sys.exit(1)

    if options.fastq is None:
        print("Error: need either single-end reads in fastq.")
        sys.exit(1)
    else:
        fastq = options.fastq

    if options.intronRef is None:
        print("Error: need intron reference file in fasta.")
        sys.exit(1)
    else:
        ref_ids, ref_seq = load_fasta(options.intronRef)
        print("reference loaded!")

    if options.out_dir is None:
        out_dir = os.path.dirname(intronRef) + "/lariatMap"
    else:
        out_dir = options.out_dir
    try:
        os.stat(out_dir)
    except:
        os.mkdir(out_dir)

    nproc = int(options.nproc)
    end_cut = int(options.end_cut)
    max_map = int(options.max_map)
    overhang = int(options.overhang)
    mismatch = int(options.mismatch)

    global FID1, FID2
    FID2 = open(out_dir + "/mapped_reads.fq", "w")
    FID1 = open(out_dir + "/lariat_reads.tab", "w")
    headline = "intron_id\tintron_order\tbp_5ss\t"
    headline += "bp_3ss\tmismatch\trlen_bp\trlen_5ss\tread_id"
    FID1.writelines(headline + "\n")
    
    cnt = 0
    if int(options.nproc) <= 1:
        aRead = []
        with open(fastq, "r") as infile:
            for line in infile:
                cnt += 1
                aRead.append(line.rstrip())
                if cnt % 4 == 0:
                    RV = map_lariat_reads(aRead, ref_ids=ref_ids, 
                        ref_seq=ref_seq, mismatch=mismatch, 
                        overhang=overhang, end_cut=end_cut, 
                        max_map=max_map)
                    aRead = []
                    show_progress(RV)
    else:
        pool = multiprocessing.Pool(processes=nproc)
        aRead = []
        with open(fastq, "r") as infile:
            for line in infile:
                cnt += 1
                aRead.append(line.rstrip())
                if cnt % 4 == 0:
                    pool.apply_async(map_lariat_reads, (aRead, ref_ids, 
                        ref_seq, mismatch, overhang, end_cut, max_map), 
                        callback=show_progress)
                    aRead = []
            pool.close()
            pool.join()
    FID1.close()
    FID2.close()
    print("")


if __name__ == "__main__":
    main()
