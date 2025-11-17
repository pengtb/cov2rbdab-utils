import requests, json
import pandas as pd

def retrieve_entity_id(pdbcode):
    baseurl = "https://data.rcsb.org/graphql?query="
    pdbcode = pdbcode.upper()
    query = "{entries(entry_ids:[\""+pdbcode+"\"]){rcsb_entry_container_identifiers{entity_ids}}}"
    url = baseurl + query
    r = requests.get(url)
    result = json.loads(r.text)
    try:
        return result['data']['entries'][0]['rcsb_entry_container_identifiers']['entity_ids']
    except:
        return result

def retrieve_instance_id(entity_id):
    baseurl = "https://data.rcsb.org/graphql?query="
    entity_id = entity_id.upper()
    query = "{polymer_entities(entity_ids:[\""+entity_id+"\"]){rcsb_polymer_entity_container_identifiers{asym_ids}}}"
    url = baseurl + query
    r = requests.get(url)
    result = json.loads(r.text)
    try:
        return result['data']['polymer_entities'][0]['rcsb_polymer_entity_container_identifiers']['asym_ids']
    except:
        return result

def retrieve_entity_annotation(entity_id):
    baseurl = "https://data.rcsb.org/graphql?query="
    entity_id = entity_id.upper()
    query = "{polymer_entities(entity_ids:[\""+entity_id+"\"]){rcsb_polymer_entity_annotation{annotation_id,name,type}}}"
    url = baseurl + query
    r = requests.get(url)
    result = json.loads(r.text)
    try:
        return result['data']['polymer_entities'][0]['rcsb_polymer_entity_annotation']
    except:
        return result

def retrieve_auth_instance_id(instance_id):
    baseurl = "https://data.rcsb.org/graphql?query="
    instance_id = instance_id.upper()
    query = "{polymer_entity_instances(instance_ids:[\""+instance_id+"\"]){rcsb_polymer_entity_instance_container_identifiers{auth_asym_id}}}"
    url = baseurl + query
    r = requests.get(url)
    result = json.loads(r.text)
    try:
        return result['data']['polymer_entity_instances'][0]['rcsb_polymer_entity_instance_container_identifiers']['auth_asym_id']
    except:
        return result
    
def retrieve_instance_annotation(instance_id):
    baseurl = "https://data.rcsb.org/graphql?query="
    instance_id = instance_id.upper()
    query = "{polymer_entity_instances(instance_ids:[\""+instance_id+"\"]){rcsb_polymer_instance_annotation{annotation_id,name,type}}}"
    url = baseurl + query
    r = requests.get(url)
    result = json.loads(r.text)
    try:
        return result['data']['polymer_entity_instances'][0]['rcsb_polymer_instance_annotation']
    except:
        return result
    
def build_idmapping(pdbcode):
    pdbcode = pdbcode.upper()
    # query
    entity_ids = retrieve_entity_id(pdbcode)
    if isinstance(entity_ids, list):
        query_entity_ids = [f'{pdbcode}_{entity_id}' for entity_id in entity_ids]
        instance_ids = [retrieve_instance_id(entity_id) for entity_id in query_entity_ids]
        query_instance_ids = [[f'{pdbcode}.{instance_id}' for instance_id in instance_subset_ids] for instance_subset_ids in instance_ids]
        auth_instance_ids = [[retrieve_auth_instance_id(instance_id) for instance_id in instance_subset_ids] for instance_subset_ids in query_instance_ids]
        entity_ids = [[entity_id] * len(instance_subset_ids) for entity_id, instance_subset_ids in zip(entity_ids,query_instance_ids)]
        # flatten
        entity_ids = [item for sublist in entity_ids for item in sublist]
        instance_ids = [item for sublist in instance_ids for item in sublist]
        auth_instance_ids = [item for sublist in auth_instance_ids for item in sublist]
        # build
        idmapping = pd.DataFrame(data={'entity_id':entity_ids,'instance_id':instance_ids,'auth_instance_id':auth_instance_ids})
        return idmapping
    else:
        print(f'{pdbcode} not found')
        return None
    
def retrieve_sequence(entity_id,  canonical=True):
    baseurl = "https://data.rcsb.org/graphql?query="
    entity_id = entity_id.upper()
    if canonical:
        key_name = 'pdbx_seq_one_letter_code_can'
    else:
        key_name = 'pdbx_seq_one_letter_code'
    query = "{polymer_entities(entity_ids:[\""+entity_id+"\"]){entity_poly{"+key_name+"}}}"
    url = baseurl + query
    r = requests.get(url)
    result = json.loads(r.text)
    try:
        return result['data']['polymer_entities'][0]['entity_poly'][key_name]
    except:
        return result
    
def retrieve_entity_description(entity_id):
    baseurl = "https://data.rcsb.org/graphql?query="
    entity_id = entity_id.upper()
    query = "{polymer_entities(entity_ids:[\""+entity_id+"\"]){rcsb_polymer_entity{pdbx_description}}}"
    url = baseurl + query
    r = requests.get(url)
    result = json.loads(r.text)
    try:
        return result['data']['polymer_entities'][0]['rcsb_polymer_entity']['pdbx_description']
    except:
        return result
    
def retrieve_entry_title(pdbcode):
    baseurl = "https://data.rcsb.org/graphql?query="
    pdbcode = pdbcode.upper()
    query = "{entries(entry_ids:[\""+pdbcode+"\"]){struct{title}}}"
    url = baseurl + query
    r = requests.get(url)
    result = json.loads(r.text)
    try:
        return result['data']['entries'][0]['struct']['title']
    except:
        return result
    
def generate_dllink(pdbcode, format='cif', zipped=False):
    pdbcode = pdbcode.upper()
    baseurl = "https://files.rcsb.org/download/"
    dlurl = baseurl + pdbcode + '.' + format
    if zipped:
        dlurl += '.gz'
    return dlurl