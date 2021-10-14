#!/bin/bash

#Author: Yang Li <yal054@ucsd.edu>
#File: run.sh
#Create Date: 2019-02-13 23:19:14
snakemake -p --rerun-incomplete -k -j 128 --cluster "qsub -l {cluster.ppn} -l {cluster.time} -N {params.jobname} -m a -q {cluster.queue} -o pbslog/{params.jobname}.pbs.out -e pbslog/{params.jobname}.pbs.err" --jobscript jobscript.pbs --jobname "{rulename}.{jobid}.pbs" --cluster-config cluster.json 2>run.log
