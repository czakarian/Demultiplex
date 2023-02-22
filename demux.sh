#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --mail-user='christinazakarian@gmail.com'
#SBATCH --mail-type=BEGIN,END,FAIL

conda activate bgmp_py39

file1="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
file2="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
file3="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
file4="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"

indexfile="/projects/bgmp/shared/2017_sequencing/indexes.txt"
outputdir="output_final"
statsfile="stats_final"

# file1="r1.fastq.gz"
# file2="r2.fastq.gz"
# file3="r3.fastq.gz"
# file4="r4.fastq.gz"

mkdir $outputdir

/usr/bin/time -v ./demux.py -i $indexfile -d $outputdir -s $statsfile -f $file1 $file2 $file3 $file4 

cd $outputdir

files=$(ls -1)
for file in $files; do gzip $file; done
