import pandas as pd
import numpy as np
from igfold import IgFoldRunner
from igfold.utils.pdb import save_PDB
import argparse, os
from tqdm import tqdm
from Bio import SeqIO
from io import StringIO
import sys

def parser_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_seqtable_path', type=str, help='path to input seqtable')
    parser.add_argument('--index_col', default=None, type=str, help='index column')
    parser.add_argument('--seq_cols', default='Hseq,Lseq', type=str, help='seq columns')
    parser.add_argument('--idx', default=None, type=int, help='index')
    parser.add_argument('--num_split', default=None, type=int, help='number of splits')
    parser.add_argument('--split_idx', default=None, type=int, help='split index')
    parser.add_argument('--refine', action='store_true', default=False, help='refine')
    parser.add_argument('--output_dir', type=str, help='path to output dir')
    parser.add_argument('--num_models', action='store', type=int, default=4, help='number of models')
    
    return parser

def standardize_seq(seq, standard_aa_str="ARNDCQEGHILKMFPSTWYV"):
    # upper
    seq = seq.upper()
    # only standard aa
    standard_aas = list(standard_aa_str)
    seq_list = list(seq)
    for idx, aa in enumerate(seq_list):
        if aa not in standard_aas:
            seq_list[idx] = 'A'
    return ''.join(seq_list)

def prepare_inputs(args, seqtable=None):
    # load table
    input_table_df = pd.read_table(args.input_seqtable_path) if seqtable is None else seqtable
    input_table_df.drop_duplicates(args.seq_cols.split(','), inplace=True)
    input_table_df.reset_index(inplace=True, drop=True)
    print('Loaded {} sequences'.format(len(input_table_df)))
    if args.index_col is None:
        input_table_df.reset_index(inplace=True)
        args.index_col = 'index'
    # subset according to idx
    if args.idx is not None:
        input_table_df = input_table_df.loc[input_table_df[args.index_col]==args.idx].copy()
    if args.num_split is not None:
        split_points = np.linspace(0, len(input_table_df), args.num_split+1, dtype=int)
        start_idx = split_points[args.split_idx]
        stop_idx = split_points[args.split_idx+1]
        print("Running on subset from {} to {}".format(start_idx, stop_idx))
        input_table_df = input_table_df.loc[start_idx:stop_idx]
    # seq cols
    Hchain_col, Lchain_col = args.seq_cols.split(',')
    # make seq dict
    seqs = []
    idxs = []
    for _, row in input_table_df.iterrows():
        seq_dict = {'H':standardize_seq(row[Hchain_col]), 'L':standardize_seq(row[Lchain_col])}
        seqs.append(seq_dict)
        idxs.append(row[args.index_col])

    return seqs, idxs

def check_model_seq(pdb_filepath, seq):
    pdb_seqrecs = list(SeqIO.parse(pdb_filepath, 'pdb-atom'))
    pdb_seq = ''.join([str(x.seq) for x in pdb_seqrecs])
    return pdb_seq == seq

def predict_structure(seqdicts, idxs, args, write2file=True, verbose=True):
    # initialize 
    runner = IgFoldRunner(num_models=args.num_models)
    # output dir
    output_dir = args.output_dir if str(args.output_dir).lower() != 'none' else './test_outputs'
    os.makedirs(output_dir, exist_ok=True)
    # run
    output_buffers = []
    iter_seqdicts = tqdm(zip(seqdicts, idxs), total=len(seqdicts)) if verbose else zip(seqdicts, idxs)
    for seqdict, idx in iter_seqdicts:
        output_filepath = os.path.join(output_dir, f'{idx}.pdb')
        wrong_output = True
        while wrong_output:
            try:
                output = runner.fold(output_filepath, sequences=seqdict, 
                                    do_refine=args.refine, do_renum=False)
            except RuntimeError:
                output = runner.fold(output_filepath, sequences=seqdict, 
                                    do_refine=False, do_renum=False)
            except ValueError:
                print(seqdict)
                sys.exit()
            # output 
            seq = "".join(seqdict.values())
            chains = list(seqdict.keys())
            delims = np.cumsum([len(s) for s in seqdict.values()]).tolist()
            if write2file:
                save_PDB(output_filepath, output.coords.squeeze(), seq, chains=chains, delim=delims, atoms=['N', 'CA', 'C', 'CB', 'O'])
                # check output pdb
                wrong_output = not check_model_seq(output_filepath, ''.join(seqdict.values()))
            else:
                pdb_string = save_PDB(None, output.coords.squeeze(), seq, chains=chains, delim=delims, atoms=['N', 'CA', 'C', 'CB', 'O'], write_pdb=False)
                output_buffer = StringIO(pdb_string)
                output_buffers.append(output_buffer)
                wrong_output = False
    if not write2file: return output_buffers
        
if __name__ == '__main__':
    args = parser_args().parse_args()
    if args.refine: 
        from igfold.refine.pyrosetta_ref import init_pyrosetta
        init_pyrosetta()
    seqdicts, idxs = prepare_inputs(args)
    predict_structure(seqdicts, idxs, args)