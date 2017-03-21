# This file is produce transcriptome sequence based on annotation file
# and genome sequence.

# Example: python genome_to_transcriptome.py -g genome.fa -a anno.gtf


import os
import sys
from optparse import OptionParser
from diceseq import load_annotation
from pyseqlib import FastaFile, fasta_write, rev_seq


def main():
    print("Welcome to Genome-to-Transcriptome!")

    # parse command line options
    parser = OptionParser()
    parser.add_option("--genome", "-g", dest="genome", default=None,
        help="The genome file in fasta.")
    parser.add_option("--anno_file", "-a", dest="anno_file", default=None,
        help="The annotation file in gtf format.")
    parser.add_option("--atype", dest="anno_type", default="GTF",
        help="The format of annotation file: GTF, GFF3 [default: %default].")
    parser.add_option("--outFile", "-o", dest="out_file", default=None, 
        help="The directory for output [default: $genome_dir/trans.fa].")
    parser.add_option("--lineLen", "-l", dest="lineLen", default="60", 
        help="The length of each line [default: %default].")

    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("use -h or --help for help on argument.")
        sys.exit(1)

    if options.anno_file == None:
        print("Error: need --anno_file for annotation.")
        sys.exit(1)
    else:
        sys.stdout.write("\rloading annotation file...")
        sys.stdout.flush()    
        anno = load_annotation(options.anno_file, options.anno_type)
        sys.stdout.write("\rloading annotation file... Done.\n")
        sys.stdout.flush()
        genes = anno["genes"]

    if options.genome is None:
        print("Error: need genome sequence in fasta.")
        sys.exit(1)
    else:
        fastaFile = FastaFile(options.genome)

    lineLen = int(options.lineLen)
    if options.out_file is None:
        fid = open(os.path.dirname(os.path.abspath(genome)) + "trans.fa", "w")
    else:
        fid = open(options.out_file, "w")

    for g in genes:
        for t in g.trans:
            chrom = t.chrom
            strand = t.strand
            id_use = g.geneID + "|" + t.tranID
            seq = []
            for e in range(t.exons.shape[0]):
                seq.append(fastaFile.get_seq(chrom, t.exons[e,0], t.exons[e,1]))
            if strand == "-" or strand == "-1": 
                seq_out = "".join(seq[::-1])
            else:
                seq_out = "".join(seq)
            fasta_write(fid, seq_out, id_use, lineLen)

    fid.close()


if __name__ == "__main__":
    main()