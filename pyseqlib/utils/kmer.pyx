## functions to count Kmers in a sequence

import numpy as np
from .fasta import get_kmer_all

def Kmer_encoder(seq, last_seq_codes=None, K=8, base_list='ACGU'):
    """get the octmer index
    """
    # code_dict = {'A': 0, 'C': 1, 'G': 2, 'U': 3}
    N_base = len(base_list) # 4 here
    code_dict = {}
    for i in range(len(base_list)):
        code_dict[base_list[i]] = i
    
    if last_seq_codes is None:
        seq_codes = [N_base**(K-1-x) * code_dict[seq[x]] for x in range(K)]
        seq = seq[K:]
    elif len(last_seq_codes) == K:        
        seq_codes = last_seq_codes
    else:
        print('Error: require len(last_seq_codes) = %d' %K)
        
    RT_vals = []
    for x in seq:
        seq_codes[:(K-1)] = [ii * N_base for ii in seq_codes[1:]]
        seq_codes[-1] = code_dict[x]
        RT_vals.append(sum(seq_codes))
        
    # return the code list and the last seq_code vector
    return RT_vals, seq_codes


def Kmer_scan(fasta_file, out_dir=None, K=8, n_gene=None, base_list='ACGU'):
    """Scan Octmers in a sequence
    """
    N_base = len(base_list)
    if n_gene is None:
        n_gene = 500000 # a default big value

    gene_list = []    
    Kmer_list = get_kmer_all(K, K, base_list)
    count_matrix = np.zeros((n_gene, N_base**K), dtype=np.int32)
    
    _gene_idx = -1
    last_seq_codes = None
    with open(fasta_file) as f:
        for line in f:
            if line.startswith('>'):
                _gene_idx += 1
                last_seq_codes = None
                gene_list.append(line.rstrip()[1:])
            else:
                seq = line.rstrip()
                _code_list, last_seq_codes = Kmer_encoder(
                    seq, last_seq_codes, K=K)
                count_matrix[_gene_idx, _code_list] += 1
                
                #seq_code_list += _code_list
                #print(_idx_list)
    # print(len(seq_code_list), len(np.unique(seq_code_list)))
    # print(gene_list)

    ## Keep found listed genes
    count_matrix = count_matrix[:(_gene_idx+1), :]

    return count_matrix, gene_list, Kmer_list
