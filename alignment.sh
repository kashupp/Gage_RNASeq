#!/bin/bash

#Path where fasta files are stored/saved.
filepath="/crex/proj/snicprojectnumber/nobackup/BEA21P160_MK"

#Path where refernce genome is stored
ref_gen_path="/crex/proj/snicprojectnumber/nobackup/labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STARgenomes/Human/GRCh38_Ensembl99_sparseD3_sjdbOverhang99"

#Path for aligned files
mkdir $filepath/alignedfiles
aligned_path="/crex/proj/snicprojectnumber/nobackup/BEA21P160_MK/alignedfiles"

: '
 unzip the files to because star does not work withcompressed files

 Loop through all files and perform fastqc
 '
for file in "$filepath"/*.fastq.gz
 do
  #Decompress the file using 4 threads
   unpigz -p 4 "$file"
 done


##Align by the STAR aligner
for file in "$filepath"/*R1*.fastq
 do
  out_name=$(echo ${file} | cut -d'R' -f 1 | cut -d'/' -f 7)
  echo $out_name

 : '
  STAR aligner is used to align reads to refernce genome.
  Path to reference genome is provided in --genomeDIr argument.
  Paired end reads are provided with --readFilesIn, first R1 and then R2.
  output file destinationand prefix to their names are provided by --outFileNamePrefix
  Aligned file is output as BAM and sorted by coordinates and is provided by
  --outSAMtype argument.
  You can also output the alignments into transcript coordinates by --quantMode. 
  '

  STAR   --runThreadN 10   --genomeDir $ref_gen_path   --readFilesIn $file "${file/R1/R2}" --outFileNamePrefix $aligned_path/$out_name  --outReadsUnmapped Fastx   --outSAMtype BAM   SortedByCoordinate      --quantMode TranscriptomeSAM   GeneCounts

  #To save the space compress the files:
   gzip "$file"
   gzip "${file/R1/R2}"
 done
