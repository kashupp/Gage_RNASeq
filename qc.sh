#!/bin/bash
#####QC of the fasta files#######
: '
 Read the number of lines in the corresponding R1 and R2 samples
 of a corressponding sample.
 Check if there are differences.
 '

#Path where fasta files are stored/saved.
filepath="/crex/proj/snic2021-22-714/nobackup/BEA21P160_MK/"

: '
 Loop on the R1 samples, extract the matching R2 and
 output the number of lines in an easy to read way.
 '

for file in "$filepath"*R1*.fastq.gz
 do
  #Print the sample identifier
  echo ${file} | cut -d'R' -f 1 | cut -d'/' -f 7

  #Output the number of lines
  wc -l < ${file} | awk '{print "Number of lines in R1 = "$1  }'
  wc -l < ${file/R1/R2} | awk '{print "Number of lines in R2 = "$1  }'
 done

##Run FASTQC and save the output

#fastqc output directory
mkdir "$filepath"fastqcOut

#Loop through all files and perform fastqc
for file in "$filepath"*
 do
  fastqc "$file" --outdir="$filepath/fastqcOut"
 done
