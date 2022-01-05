#!/bin/bash
#####Red trimming of  files#######'

#Path where fasta files are stored/saved.
filepath="/crex/proj/snic2021-22-714/nobackup/BEA21P160_MK/"

##Run Trimmomatic and save the output

#Create output directories to save the output of trimmomatic
mkdir "$filepath"trimmedFiles
mkdir "$filepath"trimmedFiles/unpaired

#Loop through all files and perform trimmin
for file in "$filepath"*R1*.fastq.gz
 do
  # Create the names if output files
   out_name_R1=$(echo ${file} | rev  | cut -d'/' -f 1  |rev)
   echo $out_name_R1

  : '
   Trimmomatic is run in paired end mode (argument PE).
   It output two types of files
   1- Trim paired, both pairs survived trimming
   They are stored in the directory "trimmedFiles".
   2- Trim unpaired, one of the pair survived trimming
   They are stored in  directory "unpaired" inside of directory "trimmedFiles"

   TruSeq3-PE adopter files are trimmed, specified by the argumment ILLUMINACLIP.
   MINLEN specified the threshold below which it will be dropped.
   In this case this threshold is 32.
  '

   trimmomatic  PE -threads 10 -phred33  $file  "${file/R1/R2}" "$filepath"trimmedFiles/"$out_name_R1"   "$filepath"trimmedFiles/unpaired/"$out_name_R1"  "$filepath"trimmedFiles/"${out_name_R1/R1/R2}" "$filepath"trimmedFiles/unpaired/"${out_name_R1/R1/R2}"  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 MINLEN:32 
 done
