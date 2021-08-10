#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1

conda activate bgmp_py39

file1="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
file2="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
file3="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
file4="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"

/usr/bin/time -v ./demux.py -i "indexes.txt" -d "output" -f $file1 $file2 $file3 $file4 
