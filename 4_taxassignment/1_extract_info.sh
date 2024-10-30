#!/bin/bash
# adapted from code found on https://www.biostars.org/p/13452/ by Frédéric Mahé
## Download NCBI's taxonomic data and GI (GenBank ID) taxonomic
## assignation.

## Variables
NCBI="ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/"
TAXDUMP="taxdump.tar.gz"
TAXID="gi_taxid_nucl.dmp.gz"
NAMES="names.dmp"
NODES="nodes.dmp"
DMP=$(echo {citations,division,gencode,merged,delnodes}.dmp)
USELESS_FILES="${TAXDUMP} ${DMP} gc.prt readme.txt"

## Download taxdump
rm -rf ${USELESS_FILES} "${NODES}" "${NAMES}"
wget "${NCBI}${TAXDUMP}" && \
    tar zxvf "${TAXDUMP}" && \
    rm -rf ${USELESS_FILES}

## Limit search space to scientific names
grep "scientific name" "${NAMES}" > "${NAMES/.dmp/_reduced.dmp}" && \
    rm -f "${NAMES}" && \
    mv "${NAMES/.dmp/_reduced.dmp}" "${NAMES}"

## Download gi_taxid_nucl
rm -f "${TAXID/.gz/}*"
wget "${NCBI}${TAXID}" && \
    gunzip "${TAXID}"

exit 0