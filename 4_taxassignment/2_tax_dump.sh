#!/bin/bash
# adapted from code found on https://www.biostars.org/p/13452/ by Frédéric Mahé

NAMES_FILE="/path/to/names.dmp"
NODES_FILE="/path/to/nodes.dmp"
TAXID="${1}"

# Read the names.dmp file into an associative array
declare -A NAMES
while IFS=$'\t' read -r -a fields; do
    TAXID=${fields[0]}
    NAME=${fields[1]}
    NAMES["$TAXID"]=$NAME
done < "$NAMES_FILE"

# Loop until you reach the root of the taxonomy (i.e., taxid = 1)
while [[ "${TAXID}" -gt 1 ]]; do
    # Obtain the parent taxa's taxid from nodes.dmp
    PARENT=$(grep --max-count=1 "^${TAXID}"$'\t' "$NODES_FILE" | cut -f3)
    
    # Build the taxonomy path
    TAXONOMY="${NAMES[$TAXID]};${TAXONOMY}"
    
    TAXID="${PARENT}"
done

echo -e "TaxID: ${1}\nTaxonomy: ${TAXONOMY}"

exit 0

