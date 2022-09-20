#!/bin/bash -l
#SBATCH -A snicprojectnumber
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 00-10:05:00
#SBATCH -J RNASeq

# Load the required modules 
module load bioinfo-tools FastQC/0.11.5 trimmomatic/0.36 star/2.7.1a subread/1.5.2 MultiQC/1.11
./qc.sh
./trimming.sh
./alignment.sh
./readcounting.sh
./multiQC.sh

