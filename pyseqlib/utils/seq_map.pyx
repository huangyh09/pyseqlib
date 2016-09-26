# Cython codes for mapping two strings

def seq_map(seq, ref, mis_match):
    RV = []
    refLen = len(ref)
    seqLen = len(seq)
    for i in range(refLen-seqLen+1):
        mis_now = 0
        seq_now = ref[i:refLen+i]
        for j in range(seqLen):
            if mis_now > mis_match: 
                break
            elif seq[j] != seq_now[j]: 
                mis_now += 1
        if mis_now <= mis_match:
            RV.append(i)
    return RV

def map_lariat_reads(read1=None, read2=None, ref_ids=None, ref_seq=None, 
    mis_match=2, overhang=30):
    mapped_bp = []
    rlen = len(read1[1])
    for i in range(len(ref_ids)):
        idx1, idx2, idx5 = [], [], []
        ref_order1, ref_order2, ref_order5 = [], [], []
        for j in range(4):
            seq = ref_seq[i][j][:overhang]
            _idx1 = seq_map(seq, read1[1], mis_match)
            for k in _idx1:
                _tmp = seq_map(read1[1][k:], ref_seq[i][j][:rlen-k+1], 
                    mis_match)
                if len(_tmp) == 1 and _tmp[0] == 0:
                    idx1.append(k)
                    ref_order1.append(j)

        for m in range(len(idx1)):
            _idx1 = idx1[m]
            if _idx1 < overhang:
                print(ref_order1[m], _idx1, rlen-_idx1, "too short.")
                continue
            print(ref_order1[m], _idx1, rlen-_idx1, "Good one!")

            for j in range(4):
                seq = read1[1][:_idx1-2] #[:_idx1]
                _idx2 = seq_map(seq, ref_seq[i][j], mis_match)

                idx2 += _idx2
                idx5 += [_idx1] * len(_idx2)
                ref_order2 += [j] * len(_idx2)
                ref_order5 += [ref_order1[m][0]] * len(_idx2)

        for m in range(len(idx2)):
            print("Lariat found!!!")
            bp = [ref_ids[i], idx2[m], ref_order2[m], ref_order5[m]]
            bp += [idx5[m], rlen-idx5[m], read1[0]]
            mapped_bp.append(bp)
    return mapped_bp

