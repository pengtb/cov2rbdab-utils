import subprocess
import pandas as pd
import numpy as np

def run_interproscan(exec_path, query_filepath, out_filepath,
                     ncpu=20, applications='SUPERFAMILY,Gene3D,CDD,SMART,Pfam',
                     formats='TSV', verbose=True):
    cmdline = f"{exec_path} -i {query_filepath} -o {out_filepath} -cpu {ncpu} -appl {applications} -f {formats}"
    cmd = cmdline.split(' ')
    if verbose:
        cmd.append('-verbose')
    print('Running {}'.format(' '.join(cmd)))
    p = subprocess.Popen(cmd, shell=False, 
                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return p

def load_result(result_filepath, **kwargs):
    columns = ['acc','md5','length','analysis','sig_acc','sig_description',
               'start','stop','score','status','date','interpro_acc','interpro_description']
    return pd.read_csv(result_filepath, sep='\t', names=columns, **kwargs)

def extract_vdomain_result(result, chain='H', dedup=True):
    filter_masks = []
    if 'SUPERFAMILY' in result.analysis.values:
        filter_masks.append((result.sig_acc=='SSF48726') & (result.start<50))
    if 'Pfam' in result.analysis.values:
        filter_masks.append(result.sig_acc=='PF07686')
    if chain == 'H':
        if 'CDD' in result.analysis.values:
            filter_masks.append(result.sig_acc=='cd04981')
        
    elif chain == 'L':
        if 'CDD' in result.analysis.values:
            filter_masks.append((result.sig_acc=='cd04980') | (result.sig_acc=='cd04984'))
    else:
        raise ValueError('chain must be H or L')
    # filter
    fileter_mask = np.stack(filter_masks, axis=0).any(0)
    vdomain_result = result.loc[fileter_mask]
    if not dedup:
        return vdomain_result
    else:
        # dedup
        sort_key = {'CDD':1, 'Pfam':2, 'SUPERFAMILY':3}
        sort_func = lambda x: sort_key[x] if x in sort_key.keys() else 100
        sorted_result = vdomain_result.sort_values('analysis', key=np.vectorize(sort_func))
        dedup_result = sorted_result.drop_duplicates('acc')
        return dedup_result