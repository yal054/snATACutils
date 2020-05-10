#!/bin/bash

#Author: Yang Li <yal054@ucsd.edu>
#File: igv.createXML.sh
#Create Date: 2019-10-04

usage() { 
cat <<EOF
Usage: igv.createXML.sh [-h] [-p <path>] [-u <URL>] [-g <genome>] [-o <output>]" 1>&2

Description:

Options:    
    -h, --help           Print help and exit
    -p, --path           local path
    -u, --url            URL path
    -g, --genome         genome mm10/hg19
    -o, --output         output xml
EOF
    exit 1
}

while getopts ":i:p:u:g:o:h" o
do
    case "$o" in
        p)
        path=${OPTARG};;
        u)
        url=${OPTARG};;
        g)
        genome=${OPTARG};;
        o) 
        output=${OPTARG};;
        *) usage;;
    esac
done
shift $((OPTIND-1))

if [ -z "${path}" ] || [ -z "${url}" ] || [ -z "${genome}" ] || [ -z "${output}" ]; then
    usage
fi

echo -e "
<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?> 
<Session genome=\"${genome}\" hasGeneTrack=\"true\" hasSequenceTrack=\"true\" version=\"8\">

<Resources>
" | sed '1d' >> ${output}

for f in `ls ${path}`;
do
echo -e "
    <Resource path=\"${url}/${f}\"/>
" >> ${output}
done;

echo -e "
</Resources>

</Session>
" >> ${output}
