# Inorder to generate sequences related to introns, and is part of intron-X

import os
import numpy as np
from pyseqlib.utils.fasta import *

def seq_maker(intron_info, fasta_file, out_dir, kmin=1, kmax=3, ss3_range=35):
    """generate the sequence feature of introns.
    1. motif sequence of 5 splice site: -4  to +7
    2. motif sequence of branch point:  -7  to +3
    3. motif sequence of 3 splice site: -16 to +4

    4. RNA sequence of whole intron
    5. RNA sequence from 5ss to bp
    6. RNA sequence from bp to 3ss

    Parameters:
    -----------
    intron_info: array like, (N, 7)
        gene_id, intron_id, chrom_id, strand, start, stop, bp
    fasta_file: string
        the file name of genome sequence in fasta
    out_dir: string
        the directory of output

    Returns:
    --------
    RV: library
        ["seq_5ss"], ["seq_BPs"], ["seq_3ss"]

    Outputs
    -------
    $out_dir/5ss_seq.fa
    $out_dir/BPs_seq.fa
    $out_dir/3ss_seq.fa
    $out_dir/intron_seq.fa
    $out_dir/5ss_BPs_seq.fa
    $out_dir/BPs_3ss_seq.fa
    $out_dir/3ss_local_seq.fa
    """

    fastaFile = FastaFile(fasta_file)
    if not os.path.exists(out_dir):
        try:
            os.makedirs(out_dir)
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise

    fid1 = open(out_dir + "/5ss_seq.fa","w")
    fid2 = open(out_dir + "/BPs_seq.fa","w")
    fid3 = open(out_dir + "/3ss_seq.fa","w")
    fid4 = open(out_dir + "/intron_seq.fa","w")
    fid5 = open(out_dir + "/5ss_BPs_seq.fa","w")
    fid6 = open(out_dir + "/BPs_3ss_seq.fa","w")
    fid7 = open(out_dir + "/3ss_local_seq.fa","w")

    kmer_lst = get_kmer_all(kmax=kmax, kmin=kmin, seqs="AUGC")
    kmer_frq = np.zeros((intron_info.shape[0], len(kmer_lst)))

    seq_5ss = []
    seq_3ss = []
    seq_BPs = []
    for i in range(intron_info.shape[0]):
        UPs = int(intron_info[i,4])
        DNs = int(intron_info[i,5])
        BPs = int(intron_info[i,6])
        chrom = str(intron_info[i,2])
        strand = str(intron_info[i,3])

        if chrom not in fastaFile.f.references:
            if chrom.count("chr") > 0:
                chrom = "".join(chrom.split("chr")[1:])
            else:
                chrom = "chr" + chrom
        if chrom not in fastaFile.f.references:
            print("Error: can't find %s in fasta file." %(intron_info[i,2]))
            sys.exit(1)

        if strand == "+":
            _seq_5ss = fastaFile.get_seq(chrom, UPs-4,  UPs+7)
            _seq_BPs = fastaFile.get_seq(chrom, BPs-7,  BPs+3)
            _seq_3ss = fastaFile.get_seq(chrom, DNs-16, DNs+4)
            
            _seq_intron  = fastaFile.get_seq(chrom, UPs, DNs)
            _seq_5ss_BPs = fastaFile.get_seq(chrom, UPs, BPs-1)
            _seq_BPs_3ss = fastaFile.get_seq(chrom, BPs+1, DNs)
            _seq_3ss_local = fastaFile.get_seq(chrom, DNs-ss3_range, 
                                               DNs+ss3_range)
        else :
            _seq_5ss = rev_seq(fastaFile.get_seq(chrom, DNs-7,  DNs+4))
            _seq_BPs = rev_seq(fastaFile.get_seq(chrom, BPs-3,  BPs+7))
            _seq_3ss = rev_seq(fastaFile.get_seq(chrom, UPs-4, UPs+16))
            
            _seq_intron  = rev_seq(fastaFile.get_seq(chrom, UPs, DNs))
            _seq_5ss_BPs = rev_seq(fastaFile.get_seq(chrom, BPs+1, DNs))
            _seq_BPs_3ss = rev_seq(fastaFile.get_seq(chrom, UPs, BPs-1))
            _seq_3ss_local = rev_seq(fastaFile.get_seq(chrom, UPs-ss3_range, 
                                                       UPs+ss3_range))
            
        seq_5ss.append(_seq_5ss)
        seq_BPs.append(_seq_BPs)
        seq_3ss.append(_seq_3ss)

        intseq = cDNA2RNA(_seq_intron)
        for j in range(len(kmer_lst)):
            kmer_frq[i,j] = get_motif(intseq, kmer_lst[j], mode="frequency")
        
        fasta_write(fid1, _seq_5ss, intron_info[i,1], length=60)
        fasta_write(fid2, _seq_BPs, intron_info[i,1], length=60)
        fasta_write(fid3, _seq_3ss, intron_info[i,1], length=60)
        fasta_write(fid4, cDNA2RNA(_seq_intron),  intron_info[i,1], length=60)
        fasta_write(fid5, cDNA2RNA(_seq_5ss_BPs), intron_info[i,1], length=60)
        fasta_write(fid6, cDNA2RNA(_seq_BPs_3ss), intron_info[i,1], length=60)
        fasta_write(fid7, cDNA2RNA(_seq_3ss_local), intron_info[i,1], length=60)

    fid1.close()
    fid2.close()
    fid3.close()
    fid4.close()
    fid5.close()
    fid6.close()
    fid7.close()

    RV = {}
    RV["seq_5ss"] = seq_5ss
    RV["seq_BPs"] = seq_BPs
    RV["seq_3ss"] = seq_3ss

    RV["kmer_lst"] = kmer_lst
    RV["kmer_frq"] = kmer_frq
    return RV

