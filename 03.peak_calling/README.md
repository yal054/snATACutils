# 03.peak calling

### 03.0 env setting
```{bash}
path0=/path_to_project
path1=/path_to_working_dir
path2script=/projects/ps-renlab/yangli/scripts/snATACutils/bin

function loadavg {
   while [ cat /proc/loadavg | awk '{print int($1)}' -gt 50 ]; do sleep 120; date; done;
}
```

you need to create one txt file for sample info, such as:
mba.dataset.list

> sample  group   techRep donor   region
> sample1  Cortex  1       1       MOp
> sample2  Cortex  1       1       MOp
> sample3  Cortex  1       2       SSp
> ...

you will also need a txt file for listing used clusters (one-column), such as:
cluster.list

> ASC
> OPC
> OGC
> MGC
> ITL23
> ITL4
> ITL5
> ITL6
> ...

### 03.1 extracting reads
```{bash}
cd $path1
mkdir bedpe qsub log; cd bedpe
path2data=/projects/ps-renlab/yangli/projects/CEMBA/00.data

for i in `sed '1d' ../mba.dataset.list | cut -f 1`;
do echo $i;
ln -s $path2data/${i}/${i}.bedpe.gz .
done;

cd $path1
cat /projects/ps-renlab/yangli/projects/CEMBA/01.cluster/mba.cerebrum.cell2anno.txt | sed -e "s/subclass/cluster/g" > mba.cerebrum.cell2anno.txt 
echo -e "
#!/bin/bash

#PBS -q hotel
#PBS -N extract_from_bedpe
#PBS -l nodes=1:ppn=6,walltime=24:00:00
#PBS -V
#PBS -M yangericlisd@gmail.com
#PBS -m a
#PBS -A ren-group
#PBS -j oe

path0=/oasis/tscc/scratch/yangli1/projects/CEMBA
path1=/oasis/tscc/scratch/yangli1/projects/CEMBA/03.peakcall
path2script=/projects/ps-renlab/yangli/scripts/snATACutils
cd \$path1

python \$path2script/snapATAC.extract_from_bedpe.py --cluster \$path1/mba.cerebrum.cell2anno.txt --indir \$path1/bedpe/ --cpu 5 --outprefix \$path1/bedpe/mba.cerebrum
" | sed '1d' > $path1/qsub/mba.cerebrum.extract_from_bedpe.qsub
qsub $path1/qsub/mba.cerebrum.extract_from_bedpe.qsub -o $path1/log/mba.cerebrum.extract_from_bedpe.log
```

alternative approach: extracting reads in sam/bam format
```{bash}
python $path2script/snapATAC.extract_from_bam.py --cluster $path1/mba.cerebrum.cell2anno.txt --indir $path0/00.data/bam/ --outprfx $path1/bam/mba.cerebrum

for i in `cat cluster.list`;
do echo $i;
samtools view -S -b bam/mba.cerebrum.cluster.${i}.sam | samtools sort - > bam/mba.cerebrum.subclass.${i}.bam
samtools index bam/mba.cerebrum.subclass.${i}.bam &
done
```

### 03.2 peak calling
```{bash}
cd $path1
mkdir tagAlign tnf5_tag peaks bw
sed -e "1d" mba.cerebrum.cell2anno.txt | cut -f 4 | sort | uniq > cluster.list
```
- peak calling for every cluster

```{bash}
for i in `cat cluster.list`;
do echo $i;
echo -e "
#!/bin/bash

#PBS -q hotel
#PBS -N ${i}.macs2
#PBS -l nodes=1:ppn=2,walltime=12:00:00
#PBS -V
#PBS -M yangericlisd@gmail.com
#PBS -m a
#PBS -A ren-group
#PBS -j oe

path0=/oasis/tscc/scratch/yangli1/projects/CEMBA
path1=/oasis/tscc/scratch/yangli1/projects/CEMBA/03.peakcall
path2script=/projects/ps-renlab/yangli/scripts/snATACutils

. ~/.bashrc
cd \$path1

zcat bedpe/mba.cerebrum.cluster.${i}.bedpe.gz | awk 'BEGIN{OFS=\"\\\t\"}{printf \"%s\\\t%s\\\t%s\\\tN\\\t1000\\\t%s\\\n%s\\\t%s\\\t%s\\\tN\\\t1000\\\t%s\\\n\",\$1,\$2,\$3,\$9,\$4,\$5,\$6,\$10}' | gzip -c > tagAlign/mba.cerebrum.cluster.${i}.tagAlign.gz
zcat tagAlign/mba.cerebrum.cluster.${i}.tagAlign.gz | awk 'BEGIN{OFS=FS=\"\\\t\"}{ if (\$6 == \"+\") {\$2 = \$2 + 4} else if (\$6 == \"-\") {\$3 = \$3 - 5} print \$0}' | gzip -nc > tnf5_tag/mba.cerebrum.cluster.${i}.tnf5_tag.gz
/home/yangli1/apps/anaconda2/bin/macs2 callpeak -t tnf5_tag/mba.cerebrum.cluster.${i}.tnf5_tag.gz -f BED -n ${i} -g hs -q 0.01 --nomodel --shift -75 --extsize 150 -B --SPMR --call-summits --outdir peaks/
grep 'chr'  peaks/${i}_treat_pileup.bdg | egrep -v \"chrM|chrUn|random\" | sort -k1,1 -k2,2n > bw/${i}_treat_pileup.srt.bdg
bedGraphToBigWig bw/${i}_treat_pileup.srt.bdg /projects/ps-renlab/yangli/genome/GRCh38.p12/GRCh38.p12.chrom.sizes  bw/${i}_treat_pileup.srt.bw
echo 'Done'
" | sed '1d' > $path1/qsub/${i}.peak.qsub
qsub $path1/qsub/${i}.peak.qsub -o $path1/log/${i}.peak.log -e $path1/log/${i}.peak.err
done
```

- peak calling for biological replicates

```{bash}
sed '1d' mba.dataset.list | awk '{print $1>>"mba.cerebrum.data.rep"$5}'

for i in `cat cluster.list`;
do echo $i;
echo -e "
#!/bin/bash

#PBS -q hotel
#PBS -N rep.${i}
#PBS -l nodes=1:ppn=1,walltime=24:00:00
#PBS -V
#PBS -M yangericlisd@gmail.com
#PBS -m a
#PBS -A ren-group
#PBS -j oe

path0=/oasis/tscc/scratch/yangli1/projects/CEMBA
path1=/oasis/tscc/scratch/yangli1/projects/CEMBA/03.peakcall
path2script=/projects/ps-renlab/yangli/scripts/snATACutils

. ~/.bashrc
cd \$path1

zcat bedpe/mba.cerebrum.cluster.${i}.bedpe.gz | grep -f mba.cerebrum.data.rep1 | gzip -c > bedpe/${i}.rep1.bedpe.gz
zcat bedpe/mba.cerebrum.cluster.${i}.bedpe.gz | grep -f mba.cerebrum.data.rep2 | gzip -c > bedpe/${i}.rep2.bedpe.gz

zcat bedpe/${i}.rep1.bedpe.gz | awk 'BEGIN{OFS=\"\\\t\"}{printf \"%s\\\t%s\\\t%s\\\tN\\\t1000\\\t%s\\\n%s\\\t%s\\\t%s\\\tN\\\t1000\\\t%s\\\n\",\$1,\$2,\$3,\$9,\$4,\$5,\$6,\$10}' | gzip -c > tagAlign/${i}.rep1.tagAlign.gz
zcat tagAlign/${i}.rep1.tagAlign.gz | awk 'BEGIN{OFS=FS=\"\\\t\"}{ if (\$6 == \"+\") {\$2 = \$2 + 4} else if (\$6 == \"-\") {\$3 = \$3 - 5} print \$0}' | gzip -nc > tnf5_tag/${i}.rep1.tnf5_tag.gz
zcat bedpe/${i}.rep2.bedpe.gz | awk 'BEGIN{OFS=\"\\\t\"}{printf \"%s\\\t%s\\\t%s\\\tN\\\t1000\\\t%s\\\n%s\\\t%s\\\t%s\\\tN\\\t1000\\\t%s\\\n\",\$1,\$2,\$3,\$9,\$4,\$5,\$6,\$10}' | gzip -c > tagAlign/${i}.rep2.tagAlign.gz
zcat tagAlign/${i}.rep2.tagAlign.gz | awk 'BEGIN{OFS=FS=\"\\\t\"}{ if (\$6 == \"+\") {\$2 = \$2 + 4} else if (\$6 == \"-\") {\$3 = \$3 - 5} print \$0}' | gzip -nc > tnf5_tag/${i}.rep2.tnf5_tag.gz

/home/yangli1/apps/anaconda2/bin/macs2 callpeak -t tnf5_tag/${i}.rep1.tnf5_tag.gz -f BED -n ${i}.rep1 -g hs -q 0.01 --nomodel --shift -75 --extsize 150 -B --SPMR --keep-dup all --call-summits --outdir peaks/
/home/yangli1/apps/anaconda2/bin/macs2 callpeak -t tnf5_tag/${i}.rep2.tnf5_tag.gz -f BED -n ${i}.rep2 -g hs -q 0.01 --nomodel --shift -75 --extsize 150 -B --SPMR --keep-dup all --call-summits --outdir peaks/
" | sed '1d' > $path1/qsub/${i}.peak2rep.qsub
qsub $path1/qsub/${i}.peak2rep.qsub -o $path1/log/${i}.peak2rep.log -e $path1/log/${i}.peak2rep.err
done
```

- peak calling for  pseudo replicates

```{bash}
for i in `cat cluster.list`;
do echo $i;
echo -e "
#!/bin/bash

#PBS -q hotel
#PBS -N pseudo.${i}
#PBS -l nodes=1:ppn=1,walltime=24:00:00
#PBS -V
#PBS -M yangericlisd@gmail.com
#PBS -m a
#PBS -A ren-group
#PBS -j oe

path0=/oasis/tscc/scratch/yangli1/projects/CEMBA
path1=/oasis/tscc/scratch/yangli1/projects/CEMBA/03.peakcall
path2script=/projects/ps-renlab/yangli/scripts/snATACutils

. ~/.bashrc
cd \$path1

# ========================
# Create pseudoReplicates
FINAL_BEDPE_FILE=\"bedpe/mba.cerebrum.cluster.${i}.bedpe.gz\"
PR_PREFIX=\"bedpe/${i}\"
PR1_TA_FILE=\"tagAlign/${i}.pseudo1.tagAlign.gz\"
PR2_TA_FILE=\"tagAlign/${i}.pseudo2.tagAlign.gz\"
TNF5_1_FILE=\"tnf5_tag/${i}.pseudo1.tnf5_tag.gz\"
TNF5_2_FILE=\"tnf5_tag/${i}.pseudo2.tnf5_tag.gz\"

# Get total number of read pairs
nlines=\$(zcat \${FINAL_BEDPE_FILE} | wc -l )
nlines=\$(( (nlines + 1) / 2 ))
# Shuffle and split BEDPE file into 2 equal parts
zcat -c \${FINAL_BEDPE_FILE} | shuf | split -d -l \${nlines} - \"\${PR_PREFIX}\" # Will produce \${PR_PREFIX}00 and \${PR_PREFIX}01
# Convert read pairs to reads into standard tagAlign file
awk 'BEGIN{OFS=\"\\\t\"}{printf \"%s\\\t%s\\\t%s\\\tN\\\t1000\\\t%s\\\n%s\\\t%s\\\t%s\\\tN\\\t1000\\\t%s\\\n\",\$1,\$2,\$3,\$9,\$4,\$5,\$6,\$10}' \"\${PR_PREFIX}00\" | gzip -c > \${PR1_TA_FILE}
rm \"\${PR_PREFIX}00\"
awk 'BEGIN{OFS=\"\\\t\"}{printf \"%s\\\t%s\\\t%s\\\tN\\\t1000\\\t%s\\\n%s\\\t%s\\\t%s\\\tN\\\t1000\\\t%s\\\n\",\$1,\$2,\$3,\$9,\$4,\$5,\$6,\$10}' \"\${PR_PREFIX}01\" | gzip -c > \${PR2_TA_FILE}
rm \"\${PR_PREFIX}01\"
# shift tag
zcat -c \${PR1_TA_FILE} | awk -F $'\\\t' 'BEGIN {OFS = FS}{ if (\$6 == \"+\") {\$2 = \$2 + 4} else if (\$6 == \"-\") {\$3 = \$3 - 5} print \$0}' | gzip -c > \${TNF5_1_FILE}
zcat -c \${PR2_TA_FILE} | awk -F $'\\\t' 'BEGIN {OFS = FS}{ if (\$6 == \"+\") {\$2 = \$2 + 4} else if (\$6 == \"-\") {\$3 = \$3 - 5} print \$0}' | gzip -c > \${TNF5_2_FILE}

# MACS2
/home/yangli1/apps/anaconda2/bin/macs2 callpeak -t \${TNF5_1_FILE} -f BED -n ${i}.pseudo1 -g hs -q 0.01 --nomodel --shift -75 --extsize 150 -B --SPMR --keep-dup all --call-summits --outdir peaks/

/home/yangli1/apps/anaconda2/bin/macs2 callpeak -t \${TNF5_2_FILE} -f BED -n ${i}.pseudo2 -g hs -q 0.01 --nomodel --shift -75 --extsize 150 -B --SPMR --keep-dup all --call-summits --outdir peaks/
echo 'Done'
" | sed '1d' > $path1/qsub/${i}.peak2pseudo.qsub
qsub $path1/qsub/${i}.peak2pseudo.qsub -o $path1/log/${i}.peak2pseudo.log -e $path1/log/${i}.peak2pseudo.err
done
```

-  check files
```{bash}
for i in `cat cluster.list`;
do echo $i;
POOLED_PEAK="peaks/${i}_peaks.narrowPeak"
REP1_PEAK="peaks/${i}.rep1_peaks.narrowPeak"
REP2_PEAK="peaks/${i}.rep2_peaks.narrowPeak"
PSEUDO1_PEAK="peaks/${i}.pseudo1_peaks.narrowPeak"
PSEUDO2_PEAK="peaks/${i}.pseudo2_peaks.narrowPeak"
BW="bw/${i}_treat_pileup.srt.bw"
REP1_BW="bw/${i}.rep1_treat_pileup.srt.bw"
REP2_BW="bw/${i}.rep2_treat_pileup.srt.bw"
PSEUDO1_BW="bw/${i}.pseudo1_treat_pileup.srt.bw"
PSEUDO2_BW="bw/${i}.pseudo2_treat_pileup.srt.bw"

[ ! -f $POOLED_PEAK ] && echo "$POOLED_PEAK does not exist" >> checkfiles.txt
[ ! -f $REP1_PEAK ] && echo "$REP1_PEAK does not exist" >> checkfiles.txt
[ ! -f $REP2_PEAK ] && echo "$REP2_PEAK does not exist" >> checkfiles.txt
[ ! -f $PSEUDO1_PEAK ] && echo "$PSEUDO1_PEAK does not exist" >> checkfiles.txt
[ ! -f $PSEUDO2_PEAK ] && echo "$PSEUDO2_PEAK does not exist" >> checkfiles.txt

[ ! -f $BW ] && echo "$BW does not exist" >> checkfiles.txt
[ ! -f $REP1_BW ] && echo "$REP1_BW does not exist" >> checkfiles.txt
[ ! -f $REP2_BW ] && echo "$REP2_BW does not exist" >> checkfiles.txt
[ ! -f $PSEUDO1_BW ] && echo "$PSEUDO1_BW does not exist" >> checkfiles.txt
[ ! -f $PSEUDO2_BW ] && echo "$PSEUDO2_BW does not exist" >> checkfiles.txt

ll $POOLED_PEAK >> checkfiles.txt
ll $REP1_PEAK >> checkfiles.txt
ll $REP2_PEAK >> checkfiles.txt
ll $PSEUDO1_PEAK >> checkfiles.txt
ll $PSEUDO2_PEAK >> checkfiles.txt
done;
```

### 03.3 parsing peaks

- navie overlap

```{bash}
mkdir parsePeak

for i in `cat cluster.list`;
do echo $i;
BLACKLIST=/projects/ps-renlab/yangli/genome/GRCh38.p12/hg38.blacklist.bed.gz
POOLED_PEAK="peaks/${i}_peaks.narrowPeak"
POOLED_SUMMIT="peaks/${i}_summits.bed"
REP1_PEAK="peaks/${i}.rep1_peaks.narrowPeak"
REP2_PEAK="peaks/${i}.rep2_peaks.narrowPeak"
PSEUDO1_PEAK="peaks/${i}.pseudo1_peaks.narrowPeak"
PSEUDO2_PEAK="peaks/${i}.pseudo2_peaks.narrowPeak"
PooledInRep1AndRep2="parsePeak/${i}.PooledInRep1AndRep2.narrowPeak.gz"
PooledInPsRep1AndPsRep2="parsePeak/${i}.PooledInPsRep1AndPsRep2.narrowPeak.gz"
naivePeakList="parsePeak/${i}.naivePeakList.narrowPeak.gz"
naiveSummitList="parsePeak/${i}.naiveSummitList.bed"

### Naive overlap
# Find pooled peaks that overlap Rep1 and Rep2 where overlap is defined as the fractional overlap wrt any one of the overlapping peak pairs >= 0.5
intersectBed -wo -a ${POOLED_PEAK} -b ${REP1_PEAK} \
| awk 'BEGIN{FS=OFS="\t"}($1~/chr*/){s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort -k1,1 -k2,2n | uniq | intersectBed -wo -a stdin -b ${REP2_PEAK} \
| awk 'BEGIN{FS=OFS="\t"}($1~/chr*/){s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort -k1,1 -k2,2n | uniq | gzip -c > ${PooledInRep1AndRep2}

# Find pooled peaks that overlap PooledPseudoRep1 and PooledPseudoRep2 where overlap is defined as the fractional overlap wrt any one of the overlapping peak pairs >= 0.5
intersectBed -wo -a ${POOLED_PEAK} -b ${PSEUDO1_PEAK} \
| awk 'BEGIN{FS=OFS="\t"}($1~/chr*/){s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort -k1,1 -k2,2n | uniq \
| intersectBed -wo -a stdin -b ${PSEUDO2_PEAK} \
| awk 'BEGIN{FS=OFS="\t"}($1~/chr*/){s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort -k1,1 -k2,2n | uniq | gzip -c > ${PooledInPsRep1AndPsRep2}

# Combine peak lists
zcat ${PooledInRep1AndRep2} ${PooledInPsRep1AndPsRep2}  | sort -k1,1 -k2,2n | uniq | awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' | grep -P 'chr[0-9XY]+(?!_)' | gzip -c  > ${naivePeakList}

# Get summit 
join -1 1 -2 4 <(zcat ${naivePeakList} | cut -f 4 | sort) <(sort -k4,4 ${POOLED_SUMMIT}) -t$'\t' \
| awk 'BEGIN{FS=OFS="\t"}{print $2,$3,$4,$1,$5}' | sort -k1,1 -k2,2n > ${naiveSummitList} 

echo -e "${i}\t${path1}/${naiveSummitList}" >> mba.cerebrum.naiveSummitList.list 
done
```

- iterative overlap peak merging

```{r}
Rscript $path2script/iterative_overlap_peak_merging.R -i mba.cerebrum.naiveSummitList.list \
        -g hg38 \
        --blacklist /projects/ps-renlab/yangli/genome/GRCh38.p12/hg38.blacklist.bed.gz \
        --chromSize /projects/ps-renlab/yangli/genome/GRCh38.p12/GRCh38.p12.chrom.sizes \
        -d parsePeak/ -o mba.cerebrum

# filter by "Score per Million" (SPM) >= 5 
sed '1d' parsePeak/mba.cerebrum.filteredNfixed.union.peakSet | awk 'BEGIN{FS=OFS="\t"}($11>=5){print $1,$2,$3,$7,$6}' | sort -k1,1 -k2,2n | uniq > mba.cerebrum.filteredNfixed.union.bed
awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,"peak"NR,$5}' mba.cerebrum.filteredNfixed.union.bed > mba.cerebrum.union.peaks.bed
```

-  assign to cell subclass

```{bash}
mkdir subclass
for i in `cat cluster.list`;
do echo $i;
intersectBed -wa -a mba.cerebrum.union.peaks.bed -b <(sed -e "1d" parsePeak/${i}.filterNfixed.peakset) > subclass/mba.cerebrum.subclass.${i}.bed
done;
```






