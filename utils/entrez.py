import pandas as pd
from tqdm import tqdm
from Bio import Entrez, SeqIO
Entrez.email = ""
Entrez.api_key = ""

def SearchGenBank(query):
    # fetch record from genbank
    handle = Entrez.esearch(db="nucleotide", term=f'"{query}"', idtype='acc', usehistory="y")
    record = Entrez.read(handle)
    handle.close()
    # convert to dataframe
    if len(record['IdList']) > 0:
        result_df = pd.DataFrame(data={'genbank': record['IdList']})
        result_df.loc[:, 'query'] = query
        result_df.loc[:, 'webenv'] = record['WebEnv']
        result_df.loc[:, 'query_key'] = record['QueryKey']
        return result_df

# upload ids to genbank history
def UploadIds(ids, batch_size=500, db='nucleotide'):
    result_df = pd.DataFrame(data={'genbank': ids})
    for startidx in tqdm(range(0, len(ids), batch_size)):
        endidx = startidx + batch_size
        if endidx > len(ids):
            endidx = len(ids)
        batch_ids = result_df.iloc[startidx:endidx, 0]
        result = Entrez.read(Entrez.epost(db=db, id=';'.join(batch_ids)))
        # save webenv and query_key
        result_df.loc[startidx:endidx, 'webenv'] = result['WebEnv']
        result_df.loc[startidx:endidx, 'query_key'] = result['QueryKey']

    return result_df

# fetch records from genbank
def FetchRecords(result_df, db='nucleotide'):
    # drop duplicate webenvs & query_keys
    result_df = result_df.drop_duplicates(subset=['webenv', 'query_key'])
    # fetch records
    records = []
    for idx, row in tqdm(result_df.iterrows(), total=len(result_df)):
        handle = Entrez.efetch(db=db, rettype="gb", retmode="text", webenv=row['webenv'], query_key=row['query_key'])
        records = records + list(SeqIO.parse(handle, "gb"))
        handle.close()
    return records

# translate mab records to protein
def TranslateMab(rec):
    # if with cds feature, use it
    feat_types = [feat.type for feat in rec.features]
    if 'CDS' in feat_types:
        cds_feat = [feat for feat in rec.features if feat.type == 'CDS'][0]
        seq = cds_feat.qualifiers['translation'][0]
    # if without cds feature, use whole sequence for translation
    else:
        seq = rec.seq.translate(to_stop=True)
    return seq