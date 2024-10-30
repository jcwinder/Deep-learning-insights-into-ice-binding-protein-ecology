# search_dictionary.py

import sys
import pickle

DICTIONARY_FILE = "/path/to/taxonomy_dictionary.pkl"

def load_dictionary_from_file(filename):
    with open(filename, 'rb') as f:
        return pickle.load(f)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python search_dictionary.py TAXID")
        sys.exit(1)

    taxid = sys.argv[1]
    dictionary = load_dictionary_from_file(DICTIONARY_FILE)

    superkingdom_taxid = dictionary.get(taxid)

    if superkingdom_taxid:
        print("TaxID:", taxid)
        print("Superkingdom TaxID:", superkingdom_taxid)
    else:
        print("TaxID:", taxid)
        print("Superkingdom TaxID: NA")
