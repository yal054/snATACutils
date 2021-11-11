# 04 peak analysis


### Unsupervised clustering using non-negative matrix factorization (NMF)

- env, data and folder setting
```{bash}
path0=/path_to_project
path1=/path_to_working_dir
path2=/path_to_nmfPmat_folder
path2script=/path_to_scripts/snATACutils/bin

function loadavg {
   while [ cat /proc/loadavg | awk '{print int($1)}' -gt 50 ]; do sleep 120; date; done;
}
```

you will need as cell-by-peak matrix after normalization
```{bash}
mkdir nmfPmat; cd nmfPmat

ln -s ../parsePmat/hba.cerebrum.parsePmat.CPMnorm.tsv .
mat="hba.cerebrum.parsePmat.CPMnorm.tsv"
prefix="hba.cerebrum.parsePmat.CPMnorm"

cut -f 1 ${mat} | sed '1d' > ${prefix}.xgi
head -n 1 ${mat} | tr '\t' '\n' | sed '1d' > ${prefix}.ygi
sed '1d' ${mat} | cut -f 2- > ${prefix}.tmp
```

convert to npz format
```{python}
import numpy as np
from scipy import sparse
from scipy.sparse import save_npz, load_npz

data = np.loadtxt('hba.cerebrum.parsePmat.CPMnorm.tmp')
data_sp = sparse.csr_matrix(data)
npz_file = "hba.cerebrum.parsePmat.CPMnorm.npz"
save_npz(npz_file, data_sp)
```

- run NMF on cell subcluster

```{bash}
cd $path2/

mkdir bin qsub log res

for r in `seq 10 27`;
do echo $r;
echo -e "
#!/bin/bash

#PBS -q hotel
#PBS -N ${r}.nmf
#PBS -l nodes=1:ppn=2,walltime=12:00:00
#PBS -V
#PBS -M yangericlisd@gmail.com
#PBS -m a
#PBS -A ren-group
#PBS -j oe

path0=/projects/ps-renlab/yangli/projects/CEMBA
path1=/projects/ps-renlab/yangli/projects/CEMBA/03.peakcall
path2=/projects/ps-renlab/yangli/projects/CEMBA/03.peakcall/nmfPmat
cd \$path2
. ~/.bashrc

python \$path2/bin/nmfATAC.cluster2peak.lite.py -i hba.cerebrum.parsePmat.CPMnorm.npz -x hba.cerebrum.parsePmat.CPMnorm.xgi -y hba.cerebrum.parsePmat.CPMnorm.ygi -r ${r} -n 10 -o res/nmfPmat.MajorType.r${r}n10 &> res/nmfPmat.MajorType.r${r}n10.log
/home/yangli1/apps/anaconda3/envs/r_env/bin/Rscript \$path2/bin/nmfATAC.plotH.R -i res/nmfPmat.MajorType.r${r}n10.H.mx -o res/nmfPmat.MajorType.r${r}n10
/home/yangli1/apps/anaconda3/envs/r_env/bin/Rscript \$path2/bin/nmfATAC.plotW.R -i res/nmfPmat.MajorType.r${r}n10.W.mx -o res/nmfPmat.MajorType.r${r}n10
python \$path2/bin/nmfATAC.stat.py -m hba.cerebrum.parsePmat.CPMnorm.npz -x hba.cerebrum.parsePmat.CPMnorm.xgi -y hba.cerebrum.parsePmat.CPMnorm.ygi --basis res/nmfPmat.MajorType.r${r}n10.W.mx --coef res/nmfPmat.MajorType.r${r}n10.H.mx -c 0.2 -o res/nmfPmat.MajorType.r${r}n10
/home/yangli1/apps/anaconda3/envs/r_env/bin/Rscript \$path2/bin/nmfATAC.statBox.R -i res/nmfPmat.MajorType.r${r}n10.statH -o res/nmfPmat.MajorType.r${r}n10 >> res/nmfPmat.MajorType.r${r}n10.sta.txt
" | sed '1d' > $path2/qsub/nmfPmat.hba.${r}.qsub
qsub $path2/qsub/nmfPmat.hba.${r}.qsub -o $path2/log/nmfPmat.hba.${r}.log
done;
```

- select stable and robust rank

```{bash}
# after all the job finished, plot measurment (sparseness)
cd res
for i in cellSparse entropy; do echo $i >> sta.lite.header; done
cut -f 1 nmfPmat.MajorType.r15n10.box.sta | sed "1 s/contributes/rank/g" > sta.box.lite.header

for r in `seq 10 27`;
do echo $r;
n=10
prefix=nmfPmat.MajorType.r${r}n${n}

# calculate cell sparseness and entropy using the statH file
Rscript $path0/bin/nmfATAC.statBox.R -i nmfPmat.MajorType.r${r}n10.statH -o nmfPmat.MajorType.r${r}n10 >> nmfPmat.MajorType.r${r}n10.sta.txt
cat ${prefix}.sta.txt | sed -e "s/^/nmfPmat\t${r}\t/g" | paste - sta.lite.header > ${prefix}.sta.tmp
sed '1d' ${prefix}.box.sta | cut -f 2 | sed -e "1i ${r}" > ${prefix}.box.contributes
sed '1d' ${prefix}.box.sta | cut -f 3 | sed -e "1i ${r}" > ${prefix}.box.sparseness
sed '1d' ${prefix}.box.sta | cut -f 4 | sed -e "1i ${r}" > ${prefix}.box.entropy
sed '1d' ${prefix}.statH | sed -e "s/^/nmfPmat\t${r}\t/g" > ${prefix}.statH.tmp
done;


cat *.sta.tmp | sed -e "1i samples\tranks\tval\tstat" | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$4,$3}' > nmfPmat.MajorType.sta.txt
rm -rf *.sta.tmp
paste sta.box.lite.header *.box.contributes > nmfPmat.MajorType.contributes.sta.txt
paste sta.box.lite.header *.box.sparseness > nmfPmat.MajorType.sparseness.sta.txt
paste sta.box.lite.header *.box.entropy > nmfPmat.MajorType.entropy.sta.txt
rm -rf *.box.contributes *.box.sparseness *.box.entropy
cat *.statH.tmp | sed -e "1i samples\tranks\txgi\tindex\tclass0\tclass1\tcontributes\tsparseness\tentropy" > nmfPmat.MajorType.statH.sta.txt
rm -rf *.statH.tmp
Rscript $path2/bin/nmfATAC.plotBox.R -i nmfPmat.MajorType.statH.sta.txt -o nmfPmat.MajorType.statH.sta
```

### Identification of differentail peaks between cell types/brain regions

- peak diff between subtypes/region
here I show differential analysi between types as example

```{bash}
mkdir betweenType

# LR test
for i in `cat CEMBA.MajorType.list`;
do echo $i;
for m in binomial negbinomial;
do echo ${m}
echo -e "
#!/bin/bash

#PBS -q hotel
#PBS -N ${i}.${m}.LRtest
#PBS -l nodes=1:ppn=4,walltime=24:00:00
#PBS -V
#PBS -M yangericlisd@gmail.com
#PBS -m a
#PBS -A ren-group
#PBS -j oe

path0=/path_to_project
path1=/path_to_working_dir/peakDiff
path2script=/path_to_scripts/snATACutils/bin

. ~/.condainit
conda activate r_env
cd \$path1
# LR test
Rscript \$path2script/snapATAC.diffpeak.LRtest.R -i RData/CEMBA.*.MajorType.${i}.sub.RData -p SubType -a SubType -m ${m} -t 1 -o betweenType/diffpeak.${i}.SubType.${m}
# JS Specificity
Rscript \$path2script/snapATAC.diffpeak.JSStest.R -i RData/CEMBA.*.MajorType.${i}.sub.RData -p SubType -a SubType -l betweenType/diffpeak.${i}.SubType.${m}.LRtest.out.tsv -s 0.001 -e 0.05 -o betweenType/diffpeak.${i}.SubType.${m}
# parse Mat
Rscript \$path2script/snapATAC.diffpeak.parseMat.R -i RData/CEMBA.*.MajorType.${i}.sub.RData -p SubType -a SubType -o betweenType/diffpeak.${i}.SubType.${m}
" | sed '1d' > $path1/qsub/diffpeak.${i}.SubType.${m}.LRtest.qsub
qsub $path1/qsub/diffpeak.${i}.SubType.${m}.LRtest.qsub -o $path1/log/diffpeak.${i}.SubType.${m}.LRtest.log
done;
done;
```





