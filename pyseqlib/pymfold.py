# This program is to run mfold with a fasta file for multiple items,
# and return the minimum free energy values.
# It is based on the mfold 3.6, which can be downloaded from here:
# http://unafold.rna.albany.edu/download/mfold-3.6.tar.gz

# Follow the instruction of INSTALL in the package.
# If you don't have the root permission, you will partially fail to install,
# but compilation will be done, you you just need to add $mfold-3.6/scripts 
# into $PATH

# demo: python pymfold.py -f in_file.fasta -o out_file.txt
# Note: the default running directory is same as the fasta data.
# But you may want to change it, as there will be many temp files.
# Then using the arguments --run_dir.

# Author: Yuanhua Huang, Y.Huang@ed.ac.uk
# Date: 2015-09-19

import os
import sys
import shutil
import subprocess
from optparse import OptionParser

import pyximport; pyximport.install()
from .utils.fasta import LightFasta

# import os
# DEVNULL = open(os.devnull, 'wb')

def get_mfold(fasta_file, out_file=None, run_dir=None, verbose=True):
    """get energy score from mfold. Make sure you have installed mfold.
    """
    fastaFile  = LightFasta(fasta_file) 
    if out_file is None:
        out_file = ".".join(fasta_file.split(".")[:-1]) + ".mfold.txt"
    if run_dir is None:
        run_dir = os.path.dirname(out_file)
        # if (os.path.exists("/tmp/")):
        #     run_dir = "/tmp/"
        # else:      
    
    RV = []
    fout = open(out_file, "w")
    fout.writelines("IDs\tlength\tenergy\n")

    cur_dir = os.getcwd()
    try:
        os.stat(run_dir + "/mfold_temp")
    except:
        os.mkdir(run_dir + "/mfold_temp")
    os.chdir(run_dir + "/mfold_temp")

    for i in range(len(fastaFile.ref)):
        fid = open("temp.fa", "w")
        fid.writelines(fastaFile.seq[i])
        fid.close()

        bashCommand = "mfold SEQ=temp.fa MAX=1 MODE=bases" 
        process = subprocess.Popen(bashCommand.split(), 
            stdout=subprocess.PIPE)

        # process = subprocess.Popen(bashCommand.split(), 
        #     stdout=subprocess.PIPE, stderr=DEVNULL, shell=False)
        
        output = process.communicate()[0]
        idx = output.index("folding energy is ")
        energy = output[idx+18:].split(" ")[0]
        if verbose is True:
            print ("The minimum free folding energy for %s is: %s." 
                %(fastaFile.ref[i], energy))
        rt_line = [fastaFile.ref[i], str(len(fastaFile.seq[i])), energy]
        fout.writelines("\t".join(rt_line) + "\n")
        
        RV.append(energy)

    shutil.rmtree("../mfold_temp", ignore_errors=True)
    os.chdir(cur_dir)
    fout.close()
    
    return RV


def main():
    print "Welcome to pymfold, a python wrap for mfold v3.6."

    #part 0. parse command line options
    parser = OptionParser()
    parser.add_option("--fasta_file", "-f", dest="fasta_file",
        help="The input fasta file for the minimum free energy by mfold.")
    parser.add_option("--out_file", "-o", dest="out_file",
        default="out_file.txt", help="The out file results.")
    parser.add_option("--run_dir", "-r", dest="run_dir", default=None,
        help="The directory for the temp files.")

    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("use -h or --help for help on argument.")
        sys.exit(1)

    run_dir = options.run_dir
    out_file = options.out_file
    fasta_file = options.fasta_file

    get_mfold(fasta_file, out_file, run_dir)


if __name__ == "__main__":
    main()