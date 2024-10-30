# generate_dictionary.py

import pickle

NAMES_FILE = "/path/to/names.dmp"
DICTIONARY_FILE = "/path/to/taxonomy_dictionary.pkl" #output path
NODES_FILE = "/path/to/nodes.dmp"

def read_names_file(names_file):
    names = {}
    with open(names_file, 'r') as f:
        for line in f:
            fields = line.strip().split("\t|\t")
            taxid, name = fields[0], fields[1]
            names[taxid] = name
    return names

def read_nodes_file(nodes_file):
    nodes = {}
    with open(nodes_file, 'r') as f:
        for line in f:
            fields = line.strip().split("\t|\t")
            taxid, parent_taxid, rank = fields[0], fields[1], fields[2]
            nodes[taxid] = {
                'parent_taxid': parent_taxid,
                'rank': rank
            }
    return nodes

def assign_superkingdom_to_dictionary(names, nodes):
    superkingdoms = {}
    for taxid in names:
        superkingdom_taxid = get_superkingdom(taxid, nodes)
        superkingdoms[taxid] = superkingdom_taxid
    return superkingdoms

def save_dictionary_to_file(dictionary, filename):
    with open(filename, 'wb') as f:
        pickle.dump(dictionary, f)

def get_superkingdom(taxid, nodes):
    while taxid != '1':
        node = nodes.get(taxid)
        if node is None:
            return None
        rank = node['rank']
        if rank == 'superkingdom':
            return taxid
        taxid = node['parent_taxid']
    return None

if __name__ == "__main__":
    names = read_names_file(NAMES_FILE)
    nodes = read_nodes_file(NODES_FILE)
    superkingdoms = assign_superkingdom_to_dictionary(names, nodes)
    save_dictionary_to_file(superkingdoms, DICTIONARY_FILE)
    print(f"Dictionary with superkingdoms saved to {DICTIONARY_FILE}")
