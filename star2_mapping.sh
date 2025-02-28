#!/bin/bash

/opt/bio/bin/STAR --genomeDir /gss/home/y.antonenkova.ext/mouse_genome/index --readFilesIn /gss/home/y.antonenkova.ext/cycloheximide/mouse_chx/samples_fastq/SRX2529040${i}.fastq  --runThreadN 12 --outFileNamePrefix /gss/home/y.antonenkova.ext/cycloheximide/mouse_chx/bams/SRX2529040${i}   --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --sjdbScore 1 --limitBAMsortRAM 60000000000
