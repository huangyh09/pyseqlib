# Cython codes for mapping two strings

def seq_map(seq, ref, mis_match, max_map=1):
    RV = []
    refLen = len(ref)
    seqLen = len(seq)
    for i in range(refLen-seqLen+1):
        mis_now = 0
        seq_now = ref[i:seqLen+i]
        for j in range(seqLen):
            if mis_now > mis_match: 
                break
            elif seq[j] != seq_now[j]: 
                mis_now += 1
        if mis_now <= mis_match:
            RV.append(i)
            mis_match = mis_now
    return RV[-max_map:]

# def map_lariat_reads(read1=None, read2=None, ref_ids=None, ref_seq=None, 
#     mis_match=2, overhang=30):
#     mapped_bp = []
#     rlen = len(read1[1])
#     for i in range(len(ref_ids)):
#         idx1, idx2, idx5 = [], [], []
#         ref_order1, ref_order2, ref_order5 = [], [], []
#         for j in range(4):
#             seq = ref_seq[i][j][1:overhang+1]#ref_seq[i][j][:overhang]
#             _idx1 = seq_map(seq, read1[1], mis_match)
#             for k in _idx1:
#                 # _tmp = seq_map(read1[1][k:], ref_seq[i][j][:rlen-k+1], 
#                 #     mis_match)
#                 _tmp = seq_map(read1[1][k:], ref_seq[i][j][1:rlen-k+2], 
#                     mis_match)
#                 if len(_tmp) == 1 and _tmp[0] == 0:
#                     idx1.append(k)
#                     ref_order1.append(j)

#         for m in range(len(idx1)):
#             _idx1 = idx1[m]
#             if _idx1 < overhang:
#                 # print(ref_order1[m], _idx1, rlen-_idx1, "too short.")
#                 continue
#             # print(ref_order1[m], _idx1, rlen-_idx1, "Good one!")

#             seq = read1[1][:_idx1-2] #[:_idx1]
#             for j in range(4):
#                 _idx2 = seq_map(seq, ref_seq[i][j], mis_match)
#                 idx2 += _idx2
#                 idx5 += [_idx1] * len(_idx2)
#                 ref_order2 += [j] * len(_idx2)
#                 ref_order5 += [ref_order1[m]] * len(_idx2)


#         for m in range(len(idx2)):
#             print("Lariat found!!!")
#             bp = [ref_ids[i], idx2[m], ref_order2[m], ref_order5[m]]
#             bp += [idx5[m], rlen-idx5[m], read1[0]]
#             mapped_bp.append(bp)
#     return mapped_bp


def map_lariat_reads(aRead=None, ref_ids=None, ref_seq=None, 
    mis_match=2, overhang=15, end_cut=1, max_map=1):
    mapped_bp = []
    rlen = len(aRead[1])
    for i in range(len(ref_ids)):
        mlen = len(ref_seq[i][0])
        for j in range(4):
            # map 5ss region onto read
            _seq = ref_seq[i][j][end_cut : overhang+end_cut]
            idx5 = seq_map(_seq, aRead[1], mis_match, max_map)

            for m in range(len(idx5)):
                # too short for bp region
                if idx5[m] < overhang: continue

                # too many errors for whole 5ss region
                _seq = aRead[1][idx5[m]:]
                _ref = ref_seq[i][j][end_cut : rlen-idx5[m]+end_cut]
                if len(seq_map(_seq, _ref, mis_match, max_map)) == 0: continue
                
                # map bp region
                _seq = aRead[1][:idx5[m]-end_cut*2]
                idxB = seq_map(_seq, ref_seq[i][j], mis_match, max_map)

                # save mapped bp info
                for b in range(len(idxB)):
                    bp_5ss = idxB[b] + len(_seq) + end_cut
                    bp_info = [ref_ids[i], j, bp_5ss, mlen-bp_5ss]
                    bp_info +=[len(_seq), rlen-len(_seq)-end_cut*2, aRead[0]]
                    mapped_bp.append(bp_info)
    RV = {}
    RV["bp"] = mapped_bp
    RV["read"] = "\n".join(aRead)
    # RV["read"] = "\n%s\n%s\n%s\n%s" %(aRead[0], aRead[1], aRead[2], aRead[3])
    return RV



