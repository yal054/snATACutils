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
```{bash}
Rscript $path2script/snapATAC.snap2cb.R -i ${snapObj}.RData -o ${prefix_name}
cbTool mtx2tsv ${prefix_name}.gmat.mtx ${prefix_name}.gene.tsv ${prefix_name}.cell.tsv exprMatrix.tsv.gz
```

Load dataset to UCSC cell browser
```bash
bash $path2script/cb.createConf.sh -i ${prefix_name} -o ./
cbBuild -o /path_to_cellbrowser/
```

### Usage of scripts

- snapATAC.leiden.py
```bash
$ python bin/snapATAC.leiden.py -h

usage: snapATAC.leiden.py [-h] [-i INPUT] [-r RESOLUTION] [-o OUTPUT]

leiden on 3-cols knn sparse matrix

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --in INPUT  input 3-col file: peak, cell, ct
  -r RESOLUTION, --resolution RESOLUTION
                        resolution from 0-1
  -o OUTPUT, --out OUTPUT
                        prefix of output file
```

- snapATAC.refineCluster.R
```bash
$ Rscript bin/snapATAC.refineCluster.R -h
usage: bin/snapATAC.refineCluster.R [-h] -i INPUT -p PARTITION -o OUTPUT

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        load preprocessed RData
  -p PARTITION, --partition PARTITION
                        partition for python-leiden
  -o OUTPUT, --output OUTPUT
                        output file prefix
```

- snapATAC.consensusLeiden.py
```bash
$ python bin/snapATAC.consensusLeiden.py -h
usage: snapATAC.consensusLeiden.py [-h] [-i INPUT] [-r RESOLUTION] [--u1 U1]
                                   [--u2 U2] [-n N] [-o OUTPUT]

leiden on 3-cols knn sparse matrix

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --in INPUT  input 3-col file: peak, cell, ct
  -r RESOLUTION, --resolution RESOLUTION
                        resolution from 0-1
  --u1 U1               left interval cutoff of CDF
  --u2 U2               right interval cutoff of CDF
  -n N                  iteration of N time
  -o OUTPUT, --out OUTPUT
                        prefix of output file
```

- snapATAC.snap2cb.R
```bash
$ Rscript bin/snapATAC.snap2cb.R -h
usage: bin/snapATAC.snap2cb.R [-h] -i INPUT -o OUTPUT

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        load preprocessed RData
  -o OUTPUT, --output OUTPUT
                        output file prefix
```


