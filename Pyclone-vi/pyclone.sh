#!/bin/bash
#SBATCH --job-name=wwy_pyclone
#SBATCH --error=%x.%j.err
#SBATCH --output=%x.%j.out
#SBATCH --partition=cpu1,cpu2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64000
cd /data/person/g5/wangwy/ULM-Project/pyclone
pyclone-vi fit -i US014-LM.tsv -o US014-LM.h5 -c 30 -d beta-binomial -r 10
pyclone-vi write-results-file -i US014-LM.h5 -o US014-LM.tsv
