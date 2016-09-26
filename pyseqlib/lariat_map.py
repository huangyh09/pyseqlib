
import os
import sys
import time
import numpy as np
import multiprocessing
from optparse import OptionParser

import pyximport; pyximport.install()
from .utils.fasta import RNA2cDNA, rev_seq
from .utils.seq_map import map_lariat_reads


# FID = None
PROCESSED = 0
MAPPED_READ = 0
START_TIME = time.time()

def show_progress(RV=None):
    global PROCESSED, MAPPED_READ, START_TIME, FID
    if len(RV) > 0:
        MAPPED_READ += 1
        for bp in RV:
            print bp
            FID.writelines("\t".join(str(x) for x in bp) + "\n")
    
    PROCESSED += 1
    if PROCESSED % 10 == 0:
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
    parser.add_option("--fastq1", "-1", dest="fastq1", default=None,
        help="The fastq file for reads mate1.")
    parser.add_option("--fastq2", "-2", dest="fastq2", default=None,
        help="The fastq file for reads mate2.")
    parser.add_option("--intronRef", "-x", dest="intronRef", default=None,
        help="The fasta file for intron sequences.")
    parser.add_option("--outDir", "-o", dest="out_dir", default=None, 
        help="The directory for output [default: same as ref].")

    parser.add_option("--nproc", "-p", dest="nproc", default=20,
        help="The number of processors to use [default: $default].")
    parser.add_option("--overHang", "-c", dest="overhang", default=20,
        help="The minimum overhang [default: $default].")
    parser.add_option("--misMatch", "-e", dest="mis_match", default=5,
        help="The maximum mis-match within overhang [default: $default].")

    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("use -h or --help for help on argument.")
        sys.exit(1)

    if options.fastq1 is None and options.fastq2 is None:
        print("Error: need either single-end or paired-end reads in fastq.")
        sys.exit(1)
    else:
        fastq1 = options.fastq1
        fastq2 = options.fastq2

    if options.intronRef is None:
        print("Error: need intron reference file in fasta.")
        sys.exit(1)
    else:
        ref_ids, ref_seq = load_fasta(options.intronRef)
        print("reference loaded!")

    if options.out_dir is None:
        out_dir = os.path.dirname(intronRef) + "/lariatMap/"
    else:
        out_dir = options.out_dir

    nproc = int(options.nproc)
    overhang = int(options.overhang)
    mis_match = int(options.mis_match)

    global FID
    FID = open(out_dir + "lariat_reads.tab", "w")
    headline = "intron_id\tintron_order\tbp_5ss\t"
    headline += "bp_3ss\tstart_bp\tstop_bp\tread_id"
    FID.writelines(headline + "\n")
    
    cnt = 0
    if int(options.nproc) <= 1:
        read1 = []
        with open(fastq1, "r") as infile:
            for line in infile:
                cnt += 1
                read1.append(line.rstrip())
                if cnt % 4 == 0:
                    RV = map_lariat_reads(read1, ref_ids=ref_ids, 
                        ref_seq=ref_seq, mis_match=mis_match, 
                        overhang=overhang)
                    read1 = []
                    show_progress(RV)
    else:
        pool = multiprocessing.Pool(processes=nproc)
        read1 = []
        with open(fastq1, "r") as infile:
            for line in infile:
                cnt += 1
                read1.append(line.rstrip())
                if cnt % 4 == 0:
                    pool.apply_async(map_lariat_reads, (read1, None, ref_ids, 
                        ref_seq, mis_match, overhang), callback=show_progress)
                    read1 = []
            pool.close()
            pool.join()
    FID.close()


if __name__ == "__main__":
    main()
