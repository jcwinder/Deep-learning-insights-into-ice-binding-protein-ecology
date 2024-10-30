#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Extract metadata from unique_accessions.csv
from Bio import Entrez
from bs4 import BeautifulSoup as bs
import subprocess
import pandas as pd
import re
import glob
from pandas import json_normalize

# Input a project accession and return sample metadata for each sample in that project, mapped to sample accession
def sample_mapping(accession, email):
    Entrez.email = email
    res = Entrez.esearch(db = "sra", term = accession, retmax = 9999)
    record = Entrez.read(res)
    print(record["Count"])
    handle = Entrez.efetch(db = "sra", id = record["IdList"], retmode = "xml") 
    readlines = handle.read()
    soup = bs(readlines, 'xml')
    experiments = soup.find_all("EXPERIMENT_PACKAGE")
    sample_attributes_mapping = {}
    for experiment in experiments:
        experiment_id = experiment.find("IDENTIFIERS").find("PRIMARY_ID").text
        sample_attributes_mapping[experiment_id] = {}
        sample_attributes = experiment.find_all("SAMPLE_ATTRIBUTE")
        for sample_attribute in sample_attributes:
            tag = sample_attribute.find("TAG").text
            value = sample_attribute.find("VALUE").text
            sample_attributes_mapping[experiment_id][tag] = value
    primary_id_mapping = {}
    for experiment in experiments:
        primary_ids = experiment.find_all("PRIMARY_ID")
        experiment_id = experiment.find("IDENTIFIERS").find("PRIMARY_ID").text
        
        for primary_id in primary_ids:
            primary_id_mapping[primary_id.text] = experiment_id
    if int(record["Count"]) == len(sample_attributes_mapping):
        print("TRUE")
        return(sample_attributes_mapping, primary_id_mapping)
    else:
        print("FALSE")

# specify file location for df column with ncbi accessions 
# make filename generalisable
csv_files = glob.glob('/path/to/metadata/folder/sra_accessions.csv')

# Initialize a set to collect unique accession numbers
unique_accessions = set()

# Read each CSV file and collect accession numbers
for file in csv_files:
    df = pd.read_csv(file)
    column_name = 'accession'
    if column_name in df.columns:
        column_values = df[column_name].tolist()
        unique_accessions.update(column_values)

# Convert the set to a list if needed
unique_accession_list = list(unique_accessions)
#Need to write a loop to generate one dictionary with all results in it
prj_accession_list = unique_accession_list

concatenated_sample_mapping = {}  # Dictionary to store the concatenated sample mappings
# Dictionary to store the mapping# Loop for Sample Mapping
concatenated_sample_mapping_ls = []
concatenated_sra_numbers = []
# Loop for Sample Mapping
for prj_accession in prj_accession_list: 
    mapped_samples = sample_mapping(prj_accession, "email@email.ac.uk")
    print(prj_accession)
    concatenated_sample_mapping = mapped_samples[0]

    # Add the 'prj_accession' to the concatenated_sample_mapping dictionary
    concatenated_sample_mapping['prj_accession'] = prj_accession
    
    samp_accessions = mapped_samples[1]
    search_term = "RR"
    filtered_dict = {key: value for key, value in samp_accessions.items() if search_term in key} 
    
    sra_numbers = list(filtered_dict.keys())
    concatenated_sra_numbers.extend(sra_numbers)
    print(concatenated_sample_mapping)
    # Store the concatenated_sample_mapping in the list
    concatenated_sample_mapping_ls.append(concatenated_sample_mapping)

df = pd.DataFrame(concatenated_sample_mapping_ls)

pd.DataFrame.head(df)


tst=concatenated_sample_mapping_ls
for item in tst:
    inner_dict=item.get(list(item.keys())[0])
    inner_dict['prj_accession']=item['prj_accession']
    inner_dict['prj_accession'] = item.pop('prj_accession', None)

tst1=[pd.DataFrame(d) for d in tst]
list_of_dataframes=tst1
merged_df = pd.concat(tst1, axis=1, sort=False)
    
merged_df.to_csv('/path/to/output/folder/filename.csv', index = True)





















