#!/bin/bash
#SBATCH --job-name=wwy_cnvkit
#SBATCH --error=%x.%j.err
#SBATCH --output=%x.%j.out
#SBATCH --partition=cpu1,cpu2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64000
cd /data/person/g5/wangwy/TJUS-WES/
cnvkit.py batch TJUS-BAM/02.Bamfile/II16-36713-16.final.bam --normal TJUS-BAM/02.Bamfile/II16-36713-14.final.bam \
-m hybrid --targets /data/person/g5/wangwy/WWY-sh/V6.bed \
--fasta /data/person/g5/wangwy/WES-Tutorial/00_ref/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
--access /data/person/g5/wangwy/WWY-sh/fixed_bed.bed  \
--output-reference reference.cnn \
--output-dir /data/person/g5/wangwy/cnvkitout/  \
--diagram --scatter
