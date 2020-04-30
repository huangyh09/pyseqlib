# This utility file is to load fasta file, 

import os
import pysam
import itertools
import numpy as np

class LightFasta:
    """This class is to load and to handle fasta file, which
    must be small size, e.g. less than 100M in plain text.
    This because it load the full data; big fasta file will
    take a lot of memenory and a long time to load. For big
    fasta files, you could use pysam.FastaFile. """
    def __init__(self, fasta_file):
        fid = open(fasta_file,"r")
        all_lines = fid.readlines()
        fid.close()
        seq, self.ref, self.seq = "", [], []
        for line in all_lines:
            line = line.split("\n")[0]
            if len(self.ref) == 0 and line[0] != ">": continue
            if line.startswith(">"):
                self.ref.append(line[1:])
                if seq == "": continue
                self.seq.append(seq)
                seq = ""
            else:
                seq = seq + line
        self.seq.append(seq)

    def get_seq(self, qref, start, stop):
        """get the sequence in a given region, the start is from 1.
        You will get stop-start+1 chars."""
        try:
            idx = self.ref.index(qref.split("chr")[-1])
        except ValueError:
            try:
                idx = self.ref.index("chr" + qref.split("chr")[-1])
            except ValueError:
                print ("No reference id as the query: %s" %qref)
                return None
        try:
            RV = self.seq[idx][start-1 : stop]
        except ValueError:
            print ("Wrong start or stop position: %d, %d" %(start, stop))
            return None
        return RV


class FastaFile:
    """This FastaFile is based on pysam.FastaFile,
    but give the same coordinate as LightFasta."""
    def __init__(self, fasta_file):
        self.f = pysam.FastaFile(fasta_file)

    def get_seq(self, qref, start, stop):
        """get the sequence in a given region, the start is from 1.
        The start and stop index may still need double check."""
        return self.f.fetch(qref, start-1, stop)


def cDNA2RNA(seq):
    """transform cDNA to RNA
    
    Parameter:
    ----------
    seq: a string, denoting sequence of cDNA

    Return:
    -------
    RV: a string, denoting according RNA
    """
    _tmp = []
    _tmp[:] = seq
    for j in range(len(_tmp)):
        if _tmp[j] == "T": _tmp[j] = "U"
    RV = "".join(_tmp)
    return RV


def RNA2cDNA(seq):
    """transform RNA to cDNA
    
    Parameter:
    ----------
    seq: a string, denoting sequence of RNA

    Return:
    -------
    RV: a string, denoting according cDNA
    """
    _tmp = []
    _tmp[:] = seq
    for j in range(len(_tmp)):
        if _tmp[j] == "U": _tmp[j] = "T"
    RV = "".join(_tmp)
    return RV


def rev_seq(seq):
    """reverse a cDNA
    
    Parameter:
    ----------
    seq: a string, denoting the original cDNA

    Return:
    -------
    RV: a string, denoting the reversed cDNA
    """
    _tmp = []
    _tmp[:] = seq
    for j in range(len(_tmp)):
        if _tmp[j] == "A": _tmp[j] = "T"
        elif _tmp[j] == "T": _tmp[j] = "A"
        elif _tmp[j] == "G": _tmp[j] = "C"
        elif _tmp[j] == "C": _tmp[j] = "G"
    RV = "".join(_tmp[::-1])
    return RV


def get_motif(seq_full, motif, mode="counts"):
    """get the counts of motif in a sequence"""
    cnt = 0
    for i in range(len(seq_full)-len(motif)+1):
        if seq_full[i:i+len(motif)] == motif:
            cnt += 1
    if mode == "counts":
        return cnt
    elif mode == "frequency":
        return cnt / (len(seq_full)-len(motif)+1.0)
    elif mode == "normalized":
        return cnt / (len(seq_full)-len(motif)+1.0) / (0.25**len(motif))
    else:
        return None


def get_kmer_all(kmax=5, kmin=1, seqs="ATGC"):
    """generate kmers"""
    RV = []
    for i in range(kmin, kmax+1):
        for _seq in itertools.product(seqs, repeat=i): 
            RV.append("".join(_seq))
    return RV


def fasta_write(fid, seq, ref, length=60):
    """write sequence into a fasta file
    
    Parameter:
    ----------
    fid: a file object, the file for writing to
    seq: a string, the sequence for writing
    ref: a string, the reference id of the sequence
    length: an int, the length of sequence each line

    Return:
    -------
    None
    """
    fid.writelines(">" + ref + "\n")
    i = -1
    for i in range(int(len(seq) / length)):
        fid.writelines(seq[i*length:(i+1)*length] + "\n")
    if (i+1)*length < len(seq):
        fid.writelines(seq[(i+1)*length:] + "\n")
    return None

def motif_score(msa, pwm_msa=None):
    """calculate motif scores

    Parameters
    ----------
    msa: list or array, element is string
    pwm_msa: list or array or None, element is string
        the msa for calculating pwm, with smooth
        if None, use msa, without smooth
        
    Return:
    -------
    score: array
        normalized motif score for each sequence
        100 means the best score according to pwm
        0 means the score for null pwm (by random)
        negetive score can happen when is poorer than null pwm
    """
    motif_len = len(msa[0])
    data = np.zeros((len(msa), motif_len), dtype="str")
    for i in range(len(msa)):
        tmp = []
        tmp[:] = msa[i].upper()
        data[i,:] = tmp
    
    if pwm_msa is None: 
        pwmS = data
        pwm_add = 0.0
    else:
        pwm_add = 0.01 # for smooth the pwm
        pwmS = np.zeros((len(pwm_msa), motif_len), dtype="str")
        for i in range(len(pwm_msa)):
            tmp = []
            tmp[:] = pwm_msa[i].upper()
            pwmS[i,:] = tmp
        
    pwm = np.zeros((4, motif_len))
    for i in range(motif_len):
        pwm[0,i] = (np.sum(pwmS[:,i]=="A")+pwm_add) / (pwmS.shape[0]+pwm_add*4)
        pwm[1,i] = (np.sum(pwmS[:,i]=="T")+pwm_add) / (pwmS.shape[0]+pwm_add*4)
        pwm[2,i] = (np.sum(pwmS[:,i]=="G")+pwm_add) / (pwmS.shape[0]+pwm_add*4)
        pwm[3,i] = (np.sum(pwmS[:,i]=="C")+pwm_add) / (pwmS.shape[0]+pwm_add*4)
        
    score = np.zeros(len(msa))
    s_max = np.sum(np.log2(pwm.max(axis=0)))
    #s_min = np.sum(np.log2(pwm.min(axis=0)))
    s_min = pwm.shape[1] * np.log2(1.0/pwm.shape[0]) #random is prefered as zero
    for i in range(data.shape[0]):
        for j in range(motif_len):
            if   data[i,j] == "A": score[i] += np.log2(pwm[0, j])
            elif data[i,j] == "T": score[i] += np.log2(pwm[1, j])
            elif data[i,j] == "G": score[i] += np.log2(pwm[2, j])
            elif data[i,j] == "C": score[i] += np.log2(pwm[3, j])
    score = (score - s_min) / (s_max-s_min) * 100
    
    return score

