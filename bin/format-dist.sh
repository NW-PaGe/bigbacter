#!/bin/bash

# format-dist.sh v1.0.0
# Author: Jared Johnson, jared.johnson@doh.wa.gov

version="v1.1.0"

#----- HELP & VERSION -----#
# help message
if [ $1 == "-h" ] || [ $1 == "--help" ] || [ $1 == "-help" ]
then
    echo -e "$0 [DIST_FILE] [TREE_FILE]" && exit 0
fi

# version
if [ $1 == "-v" ] || [ $1 == "--version" ] || [ $1 == "-version" ]
then
    echo -e ${version} && exit 0
fi

#---- INPUTS ----#
DIST_FILE=$1
TREE_FILE=$2
SOURCE=$3

#---- DETERMINE COLUMNS ----#
if [[ "$SOURCE" == 'pp' ]]
then
    echo "Using PopPUNK format."
    # poppunk
    DIST_COLS='1,2,4'
else
    echo "Using BigBacter format."
    # bb-cluster
    DIST_COLS='1,2,3'
fi

#---- EXTRACT TIPS FROM TREE ----#
cat ${TREE_FILE} | tr ',' '\n' | sed 's/:.*//g' | tr -d '(\t ' > tips.txt
# echo -e "Tree samples: \n$(cat tips.txt)"

#---- SUBSET DIST FILE ----#
# create header
echo -e "id1\tid2\tdist" > dist.formatted.txt
zcat ${DIST_FILE} \
    | tr -d '"' \
    | tr ',' '\t' \
    | cut -f ${DIST_COLS} \
    | awk -v OFS='\t' 'NR==FNR{a[$1]; next} ($1 in a) && ($2 in a)' tips.txt - >> dist.formatted.txt

#---- CLEAN UP ----#
rm tips.txt

