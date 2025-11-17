import numpy as np
from esm.inverse_folding.util import filter_backbone, get_chains
from biotite.structure.io import pdbx, pdb
import biotite.database.rcsb as rcsb
from tempfile import gettempdir
import gemmi
from scipy.spatial.distance import cdist

def get_ca_cras(pdb_file):
    """get cra of CA atoms"""
    structure = gemmi.read_structure(pdb_file)
    model = structure[0]
    cras = [cra for cra in model.all() if cra.atom.name == 'CA']
    return cras

def get_cras_pos(cras):
    """get position of cras"""
    pos = np.array([cra.atom.pos.tolist() for cra in cras])
    return pos

def get_cras_seq(cras):
    """get sequence of cras' residues"""
    residue_names = [cra.residue.name for cra in cras]
    return gemmi.one_letter_code(residue_names)

def get_cras_resid(cras, labeltype='auth'):
    """get residue ids of cras' residues"""
    if labeltype == 'auth':
        residue_ids = [cra.residue.seqid.num for cra in cras]
    elif labeltype == 'label':
        residue_ids = [cra.residue.label_seq for cra in cras]
    else:
        raise ValueError('chaintype must be auth or label')
    return np.array(residue_ids)

def get_cras_chain(cras, labeltype='auth'):
    """get chain of cras' residues"""
    if labeltype == 'auth':
        chains = [cra.chain.name for cra in cras]
    elif labeltype == 'label':
        chains = [cra.residue.subchain for cra in cras]
    else:
        raise ValueError('chaintype must be auth or label')
    return chains

def calc_distance_matrix(pos):
    """calculate distance matrix"""
    return np.sqrt(np.sum((pos[:, None, :] - pos[None, :, :])**2, axis=-1))


def load_structure(fpath, chain=None, use_author_chain=False, backbone=True, extra_fields=[]):
    """
    Args:
        fpath: filepath to either pdb or cif file
        chain: the chain id or list of chain ids to load
    Returns:
        biotite.structure.AtomArray
    """
    if fpath.endswith('cif'):
        with open(fpath) as fin:
            pdbxf = pdbx.PDBxFile.read(fin)
        structure = pdbx.get_structure(pdbxf, model=1, use_author_fields=use_author_chain,
                                       extra_fields=extra_fields)
    elif fpath.endswith('pdb'):
        with open(fpath) as fin:
            pdbf = pdb.PDBFile.read(fin)
        structure = pdb.get_structure(pdbf, model=1,
                                      extra_fields=extra_fields)
    bbmask = filter_backbone(structure) if backbone else np.ones(len(structure), dtype=bool)
    structure = structure[bbmask]
    all_chains = get_chains(structure)
    if len(all_chains) == 0:
        raise ValueError('No chains found in the input file.')
    if chain is None:
        chain_ids = all_chains
    elif isinstance(chain, list):
        chain_ids = chain
    else:
        chain_ids = [chain] 
    for chain in chain_ids:
        if chain not in all_chains:
            raise ValueError(f'Chain {chain} not found in input file')
    chain_filter = [a.chain_id in chain_ids for a in structure]
    structure = structure[chain_filter]
    return structure

def fetch_structure(pdbcode, load=True, format='mmcif', **kwargs):
    pdb_filepath = rcsb.fetch(pdbcode, format, gettempdir())
    if not load:
        return pdb_filepath
    else:
        return load_structure(pdb_filepath, **kwargs)
    
def detect_rbd_contacts_matrix(ab_instance_ids, rbd_instance_id, threshold=8):
    # extract chain ids
    instance_ids = ab_instance_ids + [rbd_instance_id]
    chain_ids = [chainid.split('.')[1] for chainid in instance_ids]
    pdbcode = ab_instance_ids[0].split('.')[0]
    # load structure
    struct = fetch_structure(pdbcode, chain=chain_ids, use_author_chain=False, extra_fields=['atom_id'])
    struct = struct[struct.atom_name=='CA']
    author_struct = fetch_structure(pdbcode, chain=None, use_author_chain=True, extra_fields=['atom_id'])
    author_struct = author_struct[author_struct.atom_name=='CA']
    # seperate as ab & rbd
    rbd_chain_id = chain_ids.pop(-1)
    rbd_struct = struct[struct.chain_id==rbd_chain_id]
    ab_struct = struct[struct.chain_id!=rbd_chain_id]
    # subset for only rbd part
    rbd_struct_atomids = rbd_struct.atom_id
    author_rbd_struct = author_struct[np.isin(author_struct.atom_id, rbd_struct_atomids)]
    author_rbd_subset_struct = author_rbd_struct[np.isin(author_rbd_struct.res_id, np.arange(319,542))]
    # coordinates
    rbd_coords = author_rbd_subset_struct.coord
    rbd_resids = author_rbd_subset_struct.res_id
    ab_coords = ab_struct.coord
    # distance
    ab_rbd_dist = cdist(ab_coords, rbd_coords)
    # contacts
    rbd_contacts_boolmask = ab_rbd_dist <= threshold
    fulllen_contacts_boolmask = np.zeros((len(ab_rbd_dist), 223), dtype=int)
    fulllen_contacts_boolmask[:, np.isin(np.arange(319,542), rbd_resids)] = rbd_contacts_boolmask
    
    return fulllen_contacts_boolmask
