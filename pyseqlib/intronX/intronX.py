# This file is designed to extract features related to introns, including
# intronL, 5ss_bpL, bp_3ssL, Delta_intron, Delta_5ss_bp, Delta_bp_3ss, 
# 5ss_motif_score, bp_motif_score, 3ss_motif_score, and 1-3 mers. WebLog will
# also be produced.

# The input is the intron_table, which contains 7 columns or 6 (if
# no branch points): gene_id, intron_id, chromasome, strand, start, stop, bp.

# If lacking branch point, its related features won't be available.


import os
import sys
import subprocess
import numpy as np
from seq_maker import seq_maker
from pyseqlib.pyRNAfold import get_RNAfold
from pyseqlib.utils.fasta import motif_score
from optparse import OptionParser, OptionGroup


def main():
    print("Welcome to Intron-X!")

    # 0. parse command line options
    parser = OptionParser()
    parser.add_option("--intron_table", "-i", dest="intron_table", default=None,
        help="The intron table file for intron information.")
    parser.add_option("--fasta", "-f", dest="fasta_file", default=None,
        help="The fasta file of genome sequence.")
    parser.add_option("--out_dir", "-o", dest="out_dir",  
        default=None, help="The directory for output [default: $fasta].")

    group = OptionGroup(parser, "Optional arguments")
    group.add_option("--kmer-range", type="int", nargs=2, dest="kmer_range",
        default=[1,3], help=("The min and max K in k-mers. [default: 1 3]"))
    group.add_option("--3ss-range", type="int", nargs=0, dest="ss3_range",
        default=35, help=("The range for 3ss structure. [default: %default]"))
    group.add_option("--no-RNAfold", action="store_true", dest="no_RNAfold", 
        default=False, help="No second strucutre for intron sequences")
    group.add_option("--no-weblogo", action="store_true", dest="no_weblogo", 
        default=False, help="No weblogo figures for motifs")
    parser.add_option_group(group)

    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("use -h or --help for help on argument.")
        sys.exit(1)
    if options.intron_table is None:
        print("Error: need --intron_table for intron information.")
        sys.exit(1)
    else:
        intron_info = np.genfromtxt(options.intron_table, delimiter="\t", 
            dtype="str", skip_header=1)
    if intron_info.shape[1] == 6:
        bp_exist = False
        intron_info = np.append(intron_info, intron_info[:,5:], axis=1)
        for i in range(intron_info.shape[0]):
            intron_info[i,6] = str((int(intron_info[i,4]) + int(intron_info[i,5]))/2)
    else:
        bp_exist = True

    if options.fasta_file is None:
        print("Error: need --fasta for fasta file of genome sequence.")
        sys.exit(1)
    else:
        fasta_file = options.fasta_file

    if options.out_dir is None:
        out_dir = os.path.dirname(fasta_file) + "/intronX"
    else:
        out_dir = options.out_dir + "/intronX"

    kmin, kmax = options.kmer_range
    ss3_range = options.ss3_range

    # 1. generate sequences for intron features
    seq_motif = seq_maker(intron_info, fasta_file, out_dir+"/seq/", 
                          kmin, kmax, ss3_range)

    # 2. calculate motif scores
    seq_score1 = motif_score(seq_motif["seq_5ss"])
    seq_score2 = motif_score(seq_motif["seq_BPs"])
    seq_score3 = motif_score(seq_motif["seq_3ss"])

    # # 3. calculate secondary structure scores
    if options.no_RNAfold is False:
        fold_score1 = get_RNAfold(out_dir+"/seq/intron_seq.fa", out_file=None)
        fold_score2 = get_RNAfold(out_dir+"/seq/5ss_BPs_seq.fa", out_file=None)
        fold_score3 = get_RNAfold(out_dir+"/seq/BPs_3ss_seq.fa", out_file=None)
        fold_score4 = get_RNAfold(out_dir+"/seq/3ss_local_seq.fa", out_file=None)
    else:
        fold_score1 = ["NA"] * (intron_info.shape[0])
        fold_score2 = ["NA"] * (intron_info.shape[0])
        fold_score3 = ["NA"] * (intron_info.shape[0])
        fold_score4 = ["NA"] * (intron_info.shape[0])

    # 4. save results
    fid = open(out_dir + "/intronX.tsv", "w")
    headline = "gene_id\tintron_id\tchrom\tstrand\tstart\tstop\tbp"
    headline += "\tintronL\t5ss_bpL\tbp_3ssL\tDeltaG_intron\tDeltaG_5ss_bp"
    headline += "\tDeltaG_bp_3ss\tDeltaG_3ss\t5ss_motif\tbp_motif\t3ss_motif"
    for _seq in seq_motif["kmer_lst"]:
        headline += "\t%s" %_seq
    fid.writelines(headline + "\n")

    for i in range(intron_info.shape[0]):
        intronL = int(intron_info[i,5]) - int(intron_info[i,4]) + 1
        ss5_bpL = int(intron_info[i,6]) - int(intron_info[i,4])
        ss3_bpL = int(intron_info[i,5]) - int(intron_info[i,6])
        if intron_info[i,3] == "-":
            ss5_bpL, ss3_bpL = ss3_bpL, ss5_bpL 

        if bp_exist is True:
            aLine = "\t".join(list(intron_info[i,:]))
            aLine += "\t%d\t%d\t%d\t" %(intronL, ss5_bpL, ss3_bpL)
            aLine += "\t".join([fold_score1[i], fold_score2[i], 
                                fold_score3[i], fold_score4[i]])
            aLine += "\t%.2f\t%.2f\t%.2f" %(seq_score1[i], seq_score2[i], 
                                            seq_score3[i])
        else:
            intron_info[i,6] = "NA"
            aLine = "\t".join(list(intron_info[i,:]))
            aLine += "\t%d\tNA\tNA\t" %(intronL)
            aLine += "\t".join([fold_score1[i], "NA", "NA", fold_score4[i]])
            aLine += "\t%.2f\tNA\t%.2f" %(seq_score1[i], seq_score3[i])
        
        for j in range(seq_motif["kmer_frq"].shape[1]):
            aLine += "\t%.2e" %seq_motif["kmer_frq"][i,j]
        fid.writelines(aLine + "\n")
    fid.close()

    if bp_exist is False:
        os.remove(out_dir+"/seq/BPs_seq.fa")
        os.remove(out_dir+"/seq/5ss_BPs_seq.fa")
        os.remove(out_dir+"/seq/BPs_3ss_seq.fa")
        if options.no_RNAfold is False:
            os.remove(out_dir+"/seq/5ss_BPs_seq.RNAfold.txt")
            os.remove(out_dir+"/seq/BPs_3ss_seq.RNAfold.txt")

    # 5. generate sequence motifs figures
    Fidx = [-4, -7, -16]
    motifs = ["/5ss_seq", "/BPs_seq", "/3ss_seq"]
    labels = ["5ss_motif", "BP_motif", "3ss_motif"]
    if options.no_weblogo:
        weblogo_idx = []
    elif bp_exist is False:
        weblogo_idx = [0, 2]
    else:
        weblogo_idx = [0,1,2]
    for i in weblogo_idx:
        label = motifs[i][1:]
        f1 = out_dir+ "/seq/" + motifs[i] 
        f2 = out_dir+motifs[i] #-s large 
        bashCommand = "weblogo -f %s.fa -o %s.pdf -F pdf -i %d --resolution 300 -t %s -c classic" %(f1, f2, Fidx[i], labels[i])
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output = process.communicate()[0]
    
    print("IntronX done!")

if __name__ == "__main__":
    main()

