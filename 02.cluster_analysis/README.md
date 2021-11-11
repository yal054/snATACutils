# 02.cluster analysis

### refine the resolution of cluster

```bash
path2script=/path_to_scripts/snATACutils/bin

mkdir refineCluster
for r in `seq 0.1 0.1 1`;
do echo ${r}
python $path2script/snapATAC.leiden.py -i hba.cortex.knn.mmtx -r ${r} -o refineCluster/hba.cortex.${r}
Rscript $path2script/snapATAC.refineCluster.R -i hba.cortex.cluster.RData -p refineCluster/hba.cortex.${r}.partition.txt -o refineCluster/hba.cortex.${r}
python $path2script/snapATAC.consensusLeiden.py -i hba.cortex.knn.mmtx -r ${r} -n 100 -o refineCluster/hba.cortex.${r}
done;

r=0.5
ln -s refineCluster/hba.cortex.${r}.refineCluster.RData hba.cortex.refineCluster.RData
ln -s refineCluster/hba.cortex.${r}.refineCluster.pdf hba.cortex.refineCluster.pdf
ln -s refineCluster/hba.cortex.${r}.refineCluster.meta.txt hba.cortex.refineCluster.meta.txt
```

### export to UCSC cell browser

Dumping snapATAC R obj to matrix
```bash
Rscript $path2script/snapATAC.snap2cb.R -i ${snapObj}.RData -o ${prefix_name}
cbTool mtx2tsv ${prefix_name}.gmat.mtx ${prefix_name}.gene.tsv ${prefix_name}.cell.tsv exprMatrix.tsv.gz
```

Load dataset to UCSC cell browser
```bash
bash $path2script/cb.createConf.sh -i ${prefix_name} -o ./
cbBuild -o /path_to_cellbrowser/
```



