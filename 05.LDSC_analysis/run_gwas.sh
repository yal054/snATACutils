#!/bin/bash

###############################################
# run GWAS (ldsc)
path0=/projects/ps-renlab/yangli/projects/CEMBA/01.joint_dat/rs1cemba/gwas/L2cluster

# liftOver and merge
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToHg19.over.chain.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToMm10.over.chain.gz

# The putative enhancer regions were mapped to the human genome (hg19) using liftOver, with a strategy similar to previous reports73. Each region was required to both uniquely map to hg19, and to uniquely map back to the original region in mm10, with the requirement that >=50% of the bases in each region were mapped back to mouse after being mapped to human.

function loadavg {
    while [ `cat /proc/loadavg | awk '{print int($1)}'` -gt 50 ]; do sleep 120; date; done;
}


mkdir raw
cd raw
# cp peaks for L2 cluster
TSCC:/oasis/tscc/scratch/yangli1/projects/CEMBA/cemba/cluster.tracks/sub_class/

ln -s ../../../ovlpPeaks/liftOver/rs1cemba.50.CREs.reciprocalToHg19.bed rs1cemba.CREs.reciprocalToHg19.bed
ln -s ../../../rs1cemba.CREs.final.bed .


mkdir data
# overlap with liftOver CREs
for i in `cat rs1cemba.L2cluster.labels.list`;
do echo $i;
intersectBed -wa -F 0.5 -a raw/rs1cemba.CREs.final.bed -b raw/${i}.q01.filteredPeakList.merged.bed > data/${i}.CREs.used.bed
join -1 1 -2 4 <(cut -f 4 data/${i}.CREs.used.bed | sort | uniq) <(sort -k4,4 raw/rs1cemba.CREs.reciprocalToHg19.bed) | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$1}' | sort -k1,1 -k2,2n | uniq > data/${i}.CREs.reciprocalToHg19.bed
done

cat data/*.CREs.used.bed | cut -f 4 | sort | uniq | wc -l
cat data/*.CREs.reciprocalToHg19.bed | cut -f 4 | sort | uniq | wc -l
cd data
ln -s ../../../ovlpPeaks/liftOver/rs1cemba.50.mm10ToHg19.map2mm10.bed rs1cemba.CREs.reciprocalToHg19.bed

# using homology peaks in default cutoff
cd $path0
mkdir pbs_log pbs_script ld_score

# run ld_score
gs=/projects/ps-renlab/yangli/genome/mm10/mm10.chrom.sizes
plink=/projects/ps-renlab/yangli/resource/1000G_EUR_Phase3_plink

function loadavg {
    while [ `cat /proc/loadavg | awk '{print int($1)}'` -gt 50 ]; do sleep 120; date; done;
}

for i in data/*.CREs.reciprocalToHg19.bed
do
    k=$(basename $i .CREs.reciprocalToHg19.bed)
    echo $k
    for j in {1..22}
    do
        ~/apps/anaconda2/bin/python ~/apps/ldsc/make_annot.py \
            --bed-file $i \
            --bimfile $plink/1000G.EUR.QC.$j.bim \
            --annot-file ld_score/${k}.$j.annot.gz &
        sleep 3
        loadavg
    done
done


## qsub_make.sh

DIR=/projects/ps-renlab/yangli/projects/CEMBA/01.joint_dat/rs1cemba/gwas/L2cluster
resources=/projects/ps-renlab/yangli/resource

for i in $(ls ld_score/*.22.annot.gz | sed 's/.22.annot.gz//g')
do
    j=$(basename $i)
    cat >pbs_script/$j.pbs <<EOF
#!/bin/bash
#PBS -q hotel
#PBS -N $j.ldscore
#PBS -l nodes=1:ppn=1
#PBS -l walltime=48:00:00
#PBS -o $DIR/pbs_log/$j.ldscore.out
#PBS -e $DIR/pbs_log/$j.ldscore.err
#PBS -V
#PBS -M yangericlisd@gmail.com
#PBS -m a
#PBS -A ren-group

export PATH="/home/yangli1/apps/anaconda2/bin:\$PATH"
cd $DIR
source activate ldsc
for i in {1..22}
do
    ~/apps/anaconda2/bin/python ~/apps/ldsc/ldsc.py\\
        --l2\\
        --bfile $resources/1000G_EUR_Phase3_plink/1000G.EUR.QC.\$i\\
        --ld-wind-cm 1\\
        --annot ld_score/$j.\$i.annot.gz\\
        --thin-annot\\
        --out ld_score/$j.\$i\\
        --print-snps $resources/1000G_Phase3_hapmap3/1000G_Phase3_hapmap3_print_snps.\$i.snp
done
source deactivate
EOF
done

for i in `ls pbs_script`;
do echo $i;
qsub pbs_script/$i 
done;

# run ldsc_test
mkdir results

ls data/*.CREs.reciprocalToHg19.bed | grep -v "rs1cemba" | sed 's/.CREs.reciprocalToHg19.bed//' | while read i
do 
    j=$(basename $i)
    printf "${j}\tld_score/${j}.,ld_score/rs1cemba.\n"
done > rs1cemba.ldcts

cts_name=rs1cemba
DIR=/projects/ps-renlab/yangli/projects/CEMBA/01.joint_dat/rs1cemba/gwas/L2cluster
sumstats=/projects/ps-renlab/yangli/resource/GWAStraits
resources=/projects/ps-renlab/yangli/resource

ls $sumstats | grep -v test | while read name
do
    i=$(basename $name .sumstats.gz)
    cat >pbs_script/$i.$cts_name.pbs <<EOF
#!/bin/bash
#PBS -q hotel
#PBS -N $i.$cts_name.ldsc
#PBS -l nodes=1:ppn=2
#PBS -l walltime=36:00:00
#PBS -o $DIR/pbs_log/$i.$cts_name.ldsc.out
#PBS -e $DIR/pbs_log/$i.$cts_name.ldsc.err
#PBS -V
#PBS -M yangericlisd@gmail.com
#PBS -m a
#PBS -A ren-group

cd $DIR
export PATH="/home/yangli1/apps/anaconda2/bin:\$PATH"
source activate ldsc
/home/yangli1/apps/anaconda2/bin/python /home/yangli1/apps/ldsc/ldsc.py\\
    --h2-cts $sumstats/$i.sumstats.gz \\
    --ref-ld-chr $resources/1000G_EUR_Phase3_baseline/baseline. \\
    --out results/${i}.$cts_name \\
    --ref-ld-chr-cts $cts_name.ldcts \\
    --w-ld-chr $resources/weights_hm3_no_hla/weights.
conda deactivate
EOF
done

for i in `ls pbs_script | grep $cts_name`;
do echo $i;
qsub pbs_script/$i -o pbs_log/$i.log
done;

# summary
for set in rs1cemba
do
    for i in results/*.${set}.cell_type_results.txt
    do
        j=$(basename $i .${set}.cell_type_results.txt)
        awk -v j=$j -v s=$set 'NR!=1{print $0"\t"j"\t"s}' $i
    done
done > res.summary.tsv

