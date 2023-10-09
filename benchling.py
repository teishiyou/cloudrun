import pandas as pd
import requests
import yaml
import re
from pathlib import Path
import sys

# read config files
with open('/home/atom/benchling/config.yaml', 'r') as y:
    config = yaml.safe_load(y)

result_path      = config['DATA']['SCREEN_RESULT']
folder_url       = config['DATA']['FOLDER_URL']
template_seq_url = config['DATA']['TEMPLATE_URL']
common_name      = config['DATA']['COMMON_NAME']


API_KEY          = config['SYSTEM']['API_KEY']
TENANT_URL       = config['SYSTEM']['TENANT_URL']
API_DNASEQ       = config['SYSTEM']['API_DNASEQ']
API_ALN          = config['SYSTEM']['API_ALN']
ALN_ALGO         = config['SYSTEM']['ALN_ALGO']

print(folder_url)
# parse folder_id
folder_id = 'lib_' + re.search('^[^-]+(?=-)', folder_url.split('/')[5])[0]

# parse template_seq_url
template_seq_id = re.search('^seq[^-]+(?=-)', template_seq_url.split('/')[6])[0]

# read screen result
def getStem(x):
    return Path(x).stem

read = pd.read_excel(result_path, sheet_name='read')
read['name'] = read['readName'].map(getStem)
print('read file shape:', read.shape)

contig = pd.read_excel(result_path, sheet_name='contig')
readme = pd.read_excel(result_path, sheet_name='readme')
description_column = pd.read_excel(result_path, sheet_name='description_column')
description_value = pd.read_excel(result_path, sheet_name='description_value')


# search via sdk
from benchling_sdk.benchling import Benchling
from benchling_sdk.auth.api_key_auth import ApiKeyAuth

benchling = Benchling(url=TENANT_URL, auth_method=ApiKeyAuth(API_KEY))

pages = benchling.dna_sequences.list(folder_id=folder_id)
print(pages)
ret = pd.DataFrame()

for page in pages:
    for seq in page:
        ## print(f"name: {seq.name} id:{seq.id}")
        ## ret = ret.append({'name': seq.name, 'id': seq.id}, ignore_index=True)
        ret = pd.concat([ret, pd.DataFrame([{'name': seq.name, 'id': seq.id}])], ignore_index=True)



print(ret.shape)
read = read.merge(ret, how='left', on='name')


# create alignment
for name, group in read.groupby('contigName'):
    print(name)

    # create sequenceId
    id_list=[]
    for i, row in group.iterrows():
        id_list.append({'sequenceId': row['id']})

    # --------------
    # create json
    # --------------
    # body = {"algorithm": "mafft",
    #         "files"    : [{"sequenceId": "seq_kSjuFZkj"}, {"sequenceId": "seq_oJnbJRc0"}],
    #         "name"     : "via api python request",
    #         "templateSequenceId": "seq_o2UVmvTN"}
    body = {'algorithm': ALN_ALGO,
            'files'    : id_list,
            'name'     : name,
            'templateSequenceId': template_seq_id}

    response = requests.post(TENANT_URL+API_ALN +':create-template-alignment', auth=(API_KEY, ""), json=body)
    print(response.json())
