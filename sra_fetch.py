#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import Entrez
from bs4 import BeautifulSoup as bs
import subprocess
import pandas as pd
import re
import glob

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


csv_files = glob.glob('') # filenames for csvs with a column containing NCBI accessions
exp_table = []

for file in csv_files:
    df = pd.read_csv(file)
    exp_table.append(df)

column_name = 'accession'
accession_number = []

for df in exp_table:
    if column_name in df.columns:
        column_values = df[column_name].tolist()
        accession_number.extend(column_values)

#exclude non-ncbi accessions
search_terms = ['PRJ', 'SRR', 'SRP', 'ERR']
pattern = "|".join(search_terms)
ncbi_accession = [item for item in accession_number if re.search(pattern, str(item), re.IGNORECASE) is not None]

#split to remove comma separated values in accession list
new_accession_list = []

for accession in ncbi_accession:
    accessions = accession.split(',')
    accessions = [acc.strip() for acc in accessions] 
    new_accession_list.extend(accessions)

print(new_accession_list)

prj_accession_list = new_accession_list

exclude_strings = []
concatenated_sample_mapping = {}  # Dictionary to store the concatenated sample mappings
concatenated_sra_numbers = []  # List to store the concatenated SRA numbers

for prj_accession in prj_accession_list: 
    mapped_samples = sample_mapping(prj_accession, "xxx@xxx.ac.uk") # add an email address
    print(prj_accession)
    concatenated_sample_mapping.update(mapped_samples[0])
    samp_accessions = mapped_samples[1]
    search_term = "RR" #get the SRA accession from the project accession
    filtered_dict = {key:value for key, value in samp_accessions.items() if search_term in key} 
    sra_numbers = list(filtered_dict.keys())
    concatenated_sra_numbers.extend(sra_numbers)

filtered_sra_numbers = [sra_number for sra_number in concatenated_sra_numbers if not any(string in sra_number for string in exclude_strings)]

# generate and write a metadata table from the dictionary


for sra_id in filtered_sra_numbers:
    print ("Currently downloading: " + sra_id)
    prefetch = "prefetch " + sra_id
    print ("The command used was: " + prefetch)
    subprocess.call(prefetch, shell=True)

# this will extract the .sra files from above into a folder named 'fastq'
for sra_id in filtered_sra_numbers:
    print("Extracting .sra files for: " + sra_id)
    move_cmd = "mv /path/to/folder/sra_files" + sra_id + "/* /path/to/folder/sra_files"
    subprocess.call(move_cmd, shell=True)
