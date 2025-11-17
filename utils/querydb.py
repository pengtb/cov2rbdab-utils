import pandas as pd

def list_all_table_names(conn):
    cursor = conn.cursor()
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
    return [row[0] for row in cursor.fetchall()]

def remove_duplicated_record(record_table, sortdict=None):
    # only qualified
    record_table = record_table.loc[record_table['notab-like']==0]
    # sort & rmdup
    seqsource_sortdict = sortdict if sortdict is not None else {'genbank':0,'patent':1,'INN':2,'sup.':3,'mutation':4,'split':5,'combination':6,'pdb':7,'CovAbDab':8}
    rmdup_record = record_table.sort_values(['source'], key=lambda x: x.apply(lambda x: seqsource_sortdict[x])).drop_duplicates(['ab_idx'])
    
    return rmdup_record

def truncate2fv(seqtable, trunct2fv_table, seq_cols=['Hseq','Lseq']):
    # make collection if sequences in seqtable
    seq_cols = sorted(seq_cols) # Hseq first, then Lseq
    current_seqs = pd.concat([seqtable[col] for col in seq_cols]).reset_index(drop=True).to_frame('seq')
    current_seqs.loc[:len(current_seqs)//2, 'chain'] = 'H'
    current_seqs.loc[len(current_seqs)//2:, 'chain'] = 'L'
    current_seqs = current_seqs.drop_duplicates('seq')
    
    # only long sequences need truncation
    current_seqs = current_seqs.loc[current_seqs.seq.apply(len)>150].reset_index(drop=True)
    
    # merge table
    merged_result = pd.merge(current_seqs, trunct2fv_table, on='seq', how='left')
    
    # seqs not in db
    notinclued_seqs = merged_result.loc[merged_result.seq_vdomain.isna(), ['seq','chain']].copy()
    
    # others continue to be truncated
    merged_result.dropna(subset=['seq_vdomain'], inplace=True)
    output_seqtable = seqtable.copy()
    for idx, row in merged_result.iterrows():
        if row['chain'] == 'H':
            output_seqtable.loc[output_seqtable.Hseq==row['seq'], 'Hseq'] = row['seq_vdomain']
        elif row['chain'] == 'L':
            output_seqtable.loc[output_seqtable.Lseq==row['seq'], 'Lseq'] = row['seq_vdomain']
    
    return output_seqtable, notinclued_seqs

def add_region_label(seqtable, region_table, seq_cols=['Hseq','Lseq']):
    # make collection if sequences in seqtable
    seq_cols = sorted(seq_cols) # Hseq first, then Lseq
    current_seqs = pd.concat([seqtable[col] for col in seq_cols]).reset_index(drop=True).to_frame('seq')
    current_seqs.loc[:len(current_seqs)//2, 'chain'] = 'H'
    current_seqs.loc[len(current_seqs)//2:, 'chain'] = 'L'
    current_seqs = current_seqs.drop_duplicates('seq')
    
    # merge table
    merged_result = pd.merge(current_seqs, region_table, on='seq', how='left')
    
    # seqs not in db
    notinclued_seqs = merged_result.loc[merged_result.region.isna(), ['seq','chain']].copy()
    
    # others continue to be truncated
    merged_result.dropna(subset=['region'], inplace=True)
    output_seqtable = seqtable.copy()
    for idx, row in merged_result.iterrows():
        if row['chain'] == 'H':
            output_seqtable.loc[output_seqtable.Hseq==row['seq'], 'Hregion'] = row['region']
        elif row['chain'] == 'L':
            output_seqtable.loc[output_seqtable.Lseq==row['seq'], 'Lregion'] = row['region']
    
    return output_seqtable, notinclued_seqs