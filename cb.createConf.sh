#!/bin/bash

#Author: Yang Li <yal054@ucsd.edu>
#File: cb.createConf.sh
#Create Date: 2019-08-09

usage() { 
cat <<EOF
Usage: cb.createConf.sh [-h] [-s <45|90>] [-p <string>]" 1>&2

Description:

Options:    
    -h, --help           Print help and exit
    -i, --inPrefix       sample/dataset lable or prefix
    -h, --hub            write cbHub or not
    -o, --outDir         output dir path
EOF
    exit 1
}

while getopts ":i:o:h" o
do
    case "$o" in
        i)
        inPrefix=$OPTARG;;
        h)
        hub=true;;
        o) 
        outDir=${OPTARG};;
        *) usage;;
    esac
done
shift $((OPTIND-1))

if [ -z "${inPrefix}" ] || [ -z "${outDir}" ]; then
    usage
fi

if [ "$hub" = true ]; then
echo -e "
# ------ conf ------ #
# summary
name = \"${inPrefix}\"
priority = 10
tags = [\"snATAC\"]
shortLabel=\"${inPrefix}\"
exprMatrix=\"${inPrefix}.exprMatrix.tsv.gz\"
geneIdType=\"symbol/gencode22\"

# data
meta=\"${inPrefix}.meta.tsv\"
coords=[
    {\"file\":\"${inPrefix}.tsne.coords.tsv\", \"shortLabel\":\"t-SNE\"},
    {\"file\":\"${inPrefix}.umap.coords.tsv\", \"shortLabel\":\"UMAP\"}
]

# optional settings
clusterField=\"cluster\"
labelField=\"cluster\"

markers=[
    {\"file\":\"${inPrefix}.diffgenes.tsv\", \"shortLabel\":\"Cluster-specific markers\"}
]

hubUrl=\"http://renlab.sdsc.edu/yangli/cbHub/${inPrefix}/hub.txt\"
#colors=\"colors.tsv\"
#acroFname = \"acronyms.tsv\"

showLabels=True
radius = 4
alpha = 0.2

quickGenesFile = \"quickGenes.csv\"

# cbHub
hubName = \"${inPrefix} Hub\"
ucscDb = \"mm10\"
unit = \"RPM (0-1 scaled)\"
" | sed '1d' > ${outDir}/cellbrowser.conf

else
echo -e "
# ------ conf ------ #
# summary
name = \"${inPrefix}\"
priority = 10
tags = [\"snATAC\"]
shortLabel=\"${inPrefix}\"
exprMatrix=\"${inPrefix}.exprMatrix.tsv.gz\"
geneIdType=\"symbol/gencode22\"

# data
meta=\"${inPrefix}.meta.tsv\"
coords=[
    {\"file\":\"${inPrefix}.tsne.coords.tsv\", \"shortLabel\":\"t-SNE\"},
    {\"file\":\"${inPrefix}.umap.coords.tsv\", \"shortLabel\":\"UMAP\"}
]

# optional settings
clusterField=\"cluster\"
labelField=\"cluster\"

markers=[
    {\"file\":\"${inPrefix}.diffgenes.tsv\", \"shortLabel\":\"Cluster-specific markers\"}
]

#hubUrl=\"http://renlab.sdsc.edu/yangli/cbHub/${inPrefix}/hub.txt\"
#colors=\"colors.tsv\"
#acroFname = \"acronyms.tsv\"

showLabels=True
radius = 4
alpha = 0.2

quickGenesFile = \"quickGenes.csv\"

# cbHub
hubName = \"${inPrefix} Hub\"
ucscDb = \"mm10\"
unit = \"RPM (0-1 scaled)\"
" | sed '1d' > ${outDir}/cellbrowser.conf

fi




