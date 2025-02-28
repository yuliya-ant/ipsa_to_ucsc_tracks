#!/bin/bash
#PBS -N mapping
#PBS -l walltime=23:59:00
#PBS -l nodes=1:ppn=12
#PBS -l mem=60gb
#PBS -d .

D=/gss/dplab/data/motoneurons-diff-chx
DIR=/gss/dplab/data/motoneurons-diff-chx/fastq/batch1
OUT_DIR=/gss/dplab/data/motoneurons-diff-chx/GRCh38/bam/batch1/
IND=/gss/dplab/genomes/GRCh38/

for f in $(ls ${DIR} | sed -e 's/_R1_001.fastq.gz//g' -e 's/_R2_001.fastq.gz//g')
do 

date

echo $f
# head -n2 ${DIR}/${f}
# --readFilesCommand zcat
/opt/bio/bin/STAR --genomeDir ${IND} --readFilesIn ${DIR}/${f}_R1_001.fastq.gz ${DIR}/${f}_R2_001.fastq.gz --readFilesCommand zcat  --runThreadN 12 --outFileNamePrefix ${OUT_DIR}/${f}_   --outSAMtype BAM SortedByCoordinate --sjdbScore 1 --limitBAMsortRAM 60000000000

done

echo 'DONE'
date
