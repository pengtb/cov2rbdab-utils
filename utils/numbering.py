import pandas as pd
import numpy as np
from tqdm import tqdm
import subprocess
from Bio.Align import PairwiseAligner, substitution_matrices
from anarci import number

def GetRegion(seq):
    standardized = standardize_seq(seq)
    numbering = GetNumbering(standardized)
    region = MarkRegion(numbering, standardized)
    return ' '.join(map(str, region))

def LoadRegionResult(filename):
    """
    Loads the result of the region annotation of AbRSA.
    """
    with open(filename, 'r') as f:
        lines = f.readlines()
    region_results = []
    count = 0
    for line in tqdm(lines):
        if line.startswith('#'):
            if line.startswith("#similarity"):
                similarity = float(line.split(' ')[1].strip())
                seq_df.loc[seq_name, 'similarity'] = similarity / 100
            else:
                continue
        elif line.startswith('>'):
            if count > 0:
                region_results.append(seq_df)
            seq_name = line[1:].strip()
            seq_df = pd.DataFrame(index=[seq_name])
            count += 1
        elif line.startswith('H') or line.startswith('L') or line.startswith('-'):
            region = line.split(':')[0].strip()[2:]
            seq = line.split(':')[1].strip()
            if region == 'EXT':
                if 'FR4' in seq_df.columns:
                    seq_df.loc[seq_name, 'FR4'] = seq_df.loc[seq_name, 'FR4'] + seq
                else:
                    seq_df.loc[seq_name, 'FR1'] = seq
            elif region == 'FR1':
                if 'FR1' in seq_df.columns:
                    seq_df.loc[seq_name, 'FR1'] = seq_df.loc[seq_name, 'FR1'] + seq
                else:
                    seq_df.loc[seq_name, 'FR1'] = seq
            else:
                seq_df.loc[seq_name, region] = seq
        else:
            continue
    # for final sequence
    region_results.append(seq_df)
    return pd.concat(region_results), lines

def RunANARCI(input_file, output_file):
    """
    Runs ANARCI on the input file.
    """
    cmd = ['ANARCI', '-i', input_file, '-o', output_file, '--csv']
    subprocess.run(cmd, env={'PATH':'/anaconda/envs/bindpredict/bin:/anaconda/condabin'})

def LoadNumbering(output_file):
    result = pd.read_csv(output_file)
    # numbering results
    numbers = np.array(result.columns[13:])
    numbers_arr = np.tile(numbers, (result.shape[0], 1))
    # mask out gaps
    seqs_arr = np.array(result.iloc[:, 13:])
    masks = seqs_arr!='-'
    masked = seqs_arr.copy()
    masked[masks] = numbers_arr[masks]
    return masked

class SeqNumbering(object):
    def __init__(self, numbering) -> None:
        self.numbering = numbering
    def __len__(self):
        return len(self.numbering)

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

def GetNumbering(seq, scheme='IMGT'):
    """Gets the numbering of the sequence."""
    # process unknown amino acids
    seq = standardize_seq(seq)
    # use ANARCI to get the numbering
    results = number(seq, scheme=scheme)[0]
    if not results:
        print(seq)
        return None
    # remove gaps
    numbering = np.array([result[0][0] for result in results])
    label = np.array([result[0][1].replace(' ', '') for result in results])
    anarci_seq = np.array([result[1] for result in results])
    valid_mask = anarci_seq!='-'
    numbering = numbering[valid_mask]
    label = label[valid_mask]
    anarci_seq = anarci_seq[valid_mask]
    # align to the original sequence
    # setup the aligner
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM90")
    aligner.internal_open_gap_score = -10
    # align
    aln = aligner.align(''.join(anarci_seq), seq)[0]
    anarci_start, anarci_stop = aln.aligned[0][0]
    seq_start, seq_stop = aln.aligned[1][0]
    # subset result
    numbering = numbering[anarci_start:anarci_stop]
    label = label[anarci_start:anarci_stop]
    numbering_result = [str(num) + l for num, l in zip(numbering, label)]
    # add gaps to unaligned positions
    numbering_result = ['-'] * seq_start + numbering_result + ['-'] * (len(seq) - seq_stop)
    return numbering_result

def MarkRegion(numbering, encode=True):
    """
    Mark regions according to IMGT numbering.
    https://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html
    """
    region = []
    for num in numbering:
        if num == '-':
            region.append('-')
        else:
            if num[-1].isalpha():
                num = num[:-1]
            num = int(num)
            if num <= 26:
                region.append('FR1')
            elif (num <= 38) & (num >= 27):
                region.append('CDR1')
            elif (num <= 55) & (num >= 39):
                region.append('FR2')
            elif (num <= 65) & (num >= 56):
                region.append('CDR2')
            elif (num <= 104) & (num >= 66):
                region.append('FR3')
            elif (num <= 117) & (num >= 105):
                region.append('CDR3')
            else:
                region.append('FR4')
    if encode:
        encode_dict = {'-':0, 'FR1':1, 'CDR1':2, 'FR2':3, 'CDR2':4, 'FR3':5, 'CDR3':6, 'FR4':7}
        region = [encode_dict[r] for r in region]
    return region

def ReindexNumberingIdxs(numbering_idxs):
    """
    Reindex the numbering indexes.
    """
    # add "A" to the end of the numbering if it not ends with a letter
    numbering_idxs = [idx + 'A' if not idx[-1].isalpha() else idx for idx in numbering_idxs]
    # fill to the beginning of the numbering with "0"
    numbering_idxs = [idx.zfill(4) for idx in numbering_idxs]
    # sort and get the reverse index
    numbering_idxs = np.array(numbering_idxs)
    sort_idxs = np.argsort(numbering_idxs)
    sorted_numbering_idxs = numbering_idxs[sort_idxs]
    # remove duplicates
    dedup_numbering_idxs = np.unique(sorted_numbering_idxs)
    # reindex the numbering with 0,1,2,3...
    new_idxs = np.arange(len(dedup_numbering_idxs))
    # reindex the original numbering
    reindex = dict(zip(dedup_numbering_idxs, new_idxs))
    reindexed_numbering_idxs = [reindex[idx] for idx in numbering_idxs]
    return reindexed_numbering_idxs
    