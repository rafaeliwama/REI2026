## loading necessary libraries
import pandas as pd
import ViromeUtils as vu
import os
import sys

### This script parses the NCBI SRA run table and kaiju results to produce a final dataframe with Linnean ranks and read counts for each taxon. ###
### Requires three input files: nodes.dmp, names.dmp, and SraRunIfo.csv ###
### This script should run in the same directory as the kaiju results files with names added. Kaiju Files should end with '.kaiju.names' ###
### library: ViromeUtils.py should be in the same directory as this script. ###
### required packages: pandas and ViromeUtils ###
### python3 KaijuToDf.py nodes.dmp names.dmp NCBI_runinfo.txt <Linnean ranks> ###


## reading NCBI taxonomy files
dics = vu.getNCBIdics(sys.argv[1], sys.argv[2])
names_dict, IDnames_dict, nodesParent_dict, nodesRank_dict = dics
NCBI_run = pd.read_csv(sys.argv[3])
names_dict, IDnames_dict, nodesParent_dict, nodesRank_dict = dics


## making starting df

NCBI_parsed = NCBI_run[['Run', 'BioSample', 'TaxID', 'size_MB', 'avgLength', 'bases', 'LibraryStrategy', 'LibrarySource', 'LibraryLayout', 'avgLength' ]]


## get Linnean ranks for each SRA run based on TaxID in NCBI_parsed

LRanks = []

for i in NCBI_parsed['TaxID']:
    rankNames = []
    for u in vu.getLRanks(i, nodesParent_dict, nodesRank_dict):
        if u not in IDnames_dict:
            rankNames.append('None')
        else:
            rankNames.append(IDnames_dict[u])
    LRanks.append(rankNames)

LRanks_df = pd.DataFrame(LRanks, columns=['species', 'genus', 'family', 'order', 'subclass', 'infraclass', 'class', 'phylum', 'kingdom', 'domain'])


NCBI_ranks_df = NCBI_parsed.join(LRanks_df)

## produces read counts for each taxon in the kaiju results
empty_list = []

for i in os.listdir():
    if i.endswith('.kaiju.names'):
        kaiju_df = pd.read_csv(i, sep='\t', names=range(0, 8))
        kaiju_df = kaiju_df[kaiju_df[0] != 'U']
        empy_list = empty_list + list(vu.get_readCounts(kaiju_df, sys.argv[4]))

kaiju_all_df = pd.DataFrame(columns=empty_list)

list_dics = []

for i in os.listdir():
    if i.endswith('.kaiju.names'):
        kaiju_df = pd.read_csv(i, sep='\t', names=range(0, 8))
        kaiju_df = kaiju_df[kaiju_df[0] != 'U']
        dict_kaiju = vu.get_readCounts(kaiju_df, sys.argv[4])
        dict_kaiju['Run'] = i.split('.')[0]
        list_dics.append(dict_kaiju)


kaiju_final = pd.DataFrame.from_records(list_dics).fillna(0)

final_virome_df = NCBI_ranks_df.merge(kaiju_final, on='Run', how='left')


final_virome_df.to_csv('final_virome_df.csv', index=False, sep='\t', header=True)