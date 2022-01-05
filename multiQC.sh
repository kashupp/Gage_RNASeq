#!/bin/bash

#File path
filepath="/crex/proj/snic2021-22-714/nobackup/BEA21P160_MK"

#FastQc file path
fastqc_reports=$filepath/fastqcOut

#Trimmed file path
trimm_files=$filepath/trimmedFiles

#Star alignment file path
star_alignments=$filepath/alignedfiles

#featurecount file path
featurecount_files=$star_alignments/countedfiles

: '
 Run multiqc with the force (-f) switch on to overwrite the old files.
 Multiple locations namely, 1= fastqc, 2= star alignment and 3=gene count files
 are provided. Multiqc produces the summary report of these process and store
 it at the destination provided by swithc -o.
 '
multiqc -f $fastqc_reports $trimm_files $star_alignments $featurecount_files -o $filepath 
'
