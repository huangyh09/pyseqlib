# Cython codes for mapping two strings

from .fasta import RNA2cDNA, rev_seq

# import re
# def seq_map(str seq, str ref, int mismatch=3, int max_map=1):
#     seq_re=re.compile('|'.join(seq[:i]+'.'+seq[i+1:] for i in range(len(seq))))
#     if len(seq_re.findall(ref)) == 0:return []
#     m = seq_re.search(ref)
#     return [m.start()]


# def seq_map(str seq, str ref, int mismatch=3, int max_map=1):
#     cnt = ref.count(seq)
#     if cnt == 0: return []
#     else: return [ref.index(seq)]


import regex
def seq_map(bytes seq, bytes ref, int mismatch=3, int max_map=1):
    # {e<=1} allows 1 substitution, insertion or deletion
    # {s<=1} allows 1 substitution 
    m = regex.findall("(%s){s<=%d}" %(seq, mismatch), ref, overlapped=True)

    RV = []
    for s in m:
        mis_cnt = 0
        #TODO: sort mis_cnt
        for i in range(len(s)):
            if s[i] != seq[i]: mis_cnt += 1
        RV.append([ref.index(s), mis_cnt])
    return RV


# #TODO: need to optimize, at least 20 times faster as str.index
# def seq_map(bytes seq, bytes ref, int mismatch=3, int max_map=1):
#     """
#     Map seq into ref, and return the start index of mapped positions.
#     This allows aligning to multiple positions, with some mismatches.

#     Parameters
#     ----------
#     seq: string
#         a pattern to find
#     ref: string
#         a reference sequence to find the pattern
#     mismatch: int
#         the maximum mismatch allowed
#     max_map: int
#         the maximum mapped positions

#     Return
#     ------
#     RV: list of int
#         the index of mapped positions; empty if none mapped
#     """
#     RV = []
#     cdef int refLen, seqLen, mis_now
#     refLen = len(ref)
#     seqLen = len(seq)
#     for i in range(refLen-seqLen+1):
#         mis_now = 0
#         for j in range(seqLen):
#             if seq[j] != ref[i+j]: 
#                 mis_now += 1
#             if mis_now > mismatch: 
#                 break
#         if mis_now <= mismatch:
#             RV.append([i, mis_now])
#             mismatch = mis_now - 0
#     return RV[-max_map:]


def map_lariat_reads(aRead=None, ref_ids=None, ref_seq=None, 
    mismatch=2, overhang=15, end_cut=1, max_map=1):
    """
    Map a lariat read onto intron. 

    Parameters
    ----------
    aRead: list of 4 strings
        a read from fastq file
    ref_ids: list of strings
        the ids of all introns
    ref_seq: list of strings
        the sequneces of all introns
    mismatch: int
        the maximum mismatch allowed
    overhang: int
        the minimum overhang required
    end_cut: int
        the number of neucleartide to cut at 5ss and bp
    max_map: int
        the maximum mapped positions

    Return
    ------
    RV["bp"]: list of int
        the index of mapped positions; empty if none mapped
    RV["read"]: list of 4 strings
        the input read
    """
    rseq = aRead[1]
    reads = [rseq, rseq[::-1], rev_seq(rseq)[::-1], rev_seq(rseq)]
    mapped_bp = []
    rlen = len(aRead[1])
    for i in range(len(ref_ids)):
        mlen = len(ref_seq[i][0])
        for j in range(4):
            # map 5ss region onto read
            _seq = ref_seq[i][0][end_cut : overhang+end_cut]
            idx5 = seq_map(_seq, reads[j], mismatch, max_map)

            for m in range(len(idx5)):
                # too short for bp region
                if idx5[m][0] < overhang: continue

                # too many errors for whole 5ss region
                _seq = reads[j][idx5[m][0]:]
                _ref = ref_seq[i][0][end_cut : rlen-idx5[m][0]+end_cut]
                if len(seq_map(_seq, _ref, mismatch, max_map)) == 0: continue
                
                # map bp region
                _seq = reads[j][:idx5[m][0]-end_cut*2]
                idxB = seq_map(_seq, ref_seq[i][0], mismatch, max_map)

                # save mapped bp info
                for b in range(len(idxB)):
                    bp_mis = idxB[b][1] + idx5[m][1]
                    bp_5ss = idxB[b][0] + len(_seq) + end_cut
                    bp_info = [ref_ids[i], j, bp_5ss, mlen-bp_5ss, bp_mis]
                    bp_info +=[len(_seq), rlen-len(_seq)-end_cut*2, aRead[0]]
                    mapped_bp.append(bp_info)
    RV = {}
    RV["bp"] = mapped_bp
    RV["read"] = "\n".join(aRead)
    return RV





# # TODO: It is not good to have reverse reference as the 5ss is fixed.
# def map_lariat_reads(aRead=None, ref_ids=None, ref_seq=None, 
#     mismatch=2, overhang=15, end_cut=1, max_map=1):
#     """
#     Map a lariat read onto intron. 

#     Parameters
#     ----------
#     aRead: list of 4 strings
#         a read from fastq file
#     ref_ids: list of strings
#         the ids of all introns
#     ref_seq: list of strings
#         the sequneces of all introns
#     mismatch: int
#         the maximum mismatch allowed
#     overhang: int
#         the minimum overhang required
#     end_cut: int
#         the number of neucleartide to cut at 5ss and bp
#     max_map: int
#         the maximum mapped positions

#     Return
#     ------
#     RV["bp"]: list of int
#         the index of mapped positions; empty if none mapped
#     RV["read"]: list of 4 strings
#         the input read
#     """
#     mapped_bp = []
#     rlen = len(aRead[1])
#     for i in range(len(ref_ids)):
#         mlen = len(ref_seq[i][0])
#         for j in range(4):
#             # map 5ss region onto read
#             _seq = ref_seq[i][j][end_cut : overhang+end_cut]
#             idx5 = seq_map(_seq, aRead[1], mismatch, max_map)
#             # print len(idx5)

#             for m in range(len(idx5)):
#                 # too short for bp region
#                 if idx5[m][0] < overhang: continue

#                 print ref_ids[i], j
#                 print "G"+_seq
#                 print ref_seq[i][j][0 : 20]
#                 print aRead[1][idx5[m][0]-1:]

#                 # too many errors for whole 5ss region
#                 _seq = aRead[1][idx5[m][0]:]
#                 _ref = ref_seq[i][j][end_cut : rlen-idx5[m][0]+end_cut]
#                 if len(seq_map(_seq, _ref, mismatch, max_map)) == 0: continue
                
#                 # map bp region
#                 _seq = aRead[1][:idx5[m][0]-end_cut*2]
#                 idxB = seq_map(_seq, ref_seq[i][j], mismatch, max_map)

#                 # save mapped bp info
#                 for b in range(len(idxB)):
#                     bp_5ss = idxB[b][0] + len(_seq) + end_cut
#                     bp_info = [ref_ids[i], j, bp_5ss, mlen-bp_5ss, idxB[b][1]]
#                     bp_info +=[len(_seq), rlen-len(_seq)-end_cut*2, aRead[0]]
#                     mapped_bp.append(bp_info)
#     RV = {}
#     RV["bp"] = mapped_bp
#     RV["read"] = "\n".join(aRead)
#     return RV

