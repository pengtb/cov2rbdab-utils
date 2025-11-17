import pandas as pd
import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices

rbd_wt_seq = "RVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNF"
rbd_wt_seqarr = np.asarray(list(rbd_wt_seq))
rbd_wt_resids = np.arange(319,542)

def load_variant_rbdseq_table(seqtable_filepath):
    variant_rbdseq = pd.read_table(seqtable_filepath)
    length_table = len(variant_rbdseq)
    variant_rbdseq.loc[length_table, 'lineage'] = 'WT'
    variant_rbdseq.loc[length_table, 'rbd_seq'] = rbd_wt_seq
    return variant_rbdseq

def get_aln_score(query_seq, target_seq):
    # aligner
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load('BLOSUM90')
    # aln
    aln = aligner.align(query_seq, target_seq)
    return aln, aln.score

def identify_lineage(variant_rbdseq_table, query_rbdseq, 
                     rbdseq_colname='rbd_seq', lineage_colname='lineage'):
    # aligner
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load('BLOSUM90')
    
    # alignment
    alns = [aligner.align(query_rbdseq, rbdseq) for rbdseq in variant_rbdseq_table[rbdseq_colname].values]
    aln_scores = np.asarray([aln.score for aln in alns])
    
    # lineage with highest score
    matched_lineages = variant_rbdseq_table[lineage_colname].values[np.argmax(aln_scores)]
    
    return matched_lineages