#!/bin/bash

#Path where aligned files are stored/saved.
filepath="/crex/proj/snic2021-22-714/nobackup/BEA21P160_MK/alignedfiles"

#Path where refernce genome is stored
ref_gen_path="/crex/proj/snic2021-22-714/nobackup/labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STARgenomes/Human/GRCh38_Ensembl99_sparseD3_sjdbOverhang99"

#Path for counted files
mkdir $filepath/countedfiles


##read count  by subread package using featurecount software
for file in "$filepath"/*Aligned.sortedByCoord.out.bam
do
out_name=$(echo ${file} | rev |  cut -d'/' -f 1 | cut -d'.' -f 4 | rev)_Count

echo $out_name

: '
-p if the input datacontain paired-end read.
The input is the paired endindicated by the presence of switch -B. 
Annotation file is provided with the argument -a. It is a gtf file.
-g specifies the attribute type to group featues. Itcould be exon or gene_id
-t is the feature type.
-o output destination and file namae.
An aligned and  sorted by coordibate file is provided at the end  of command.  
'

featureCounts -T 10 -p -B -a $ref_gen_path/Homo_sapiens.GRCh38.99.gtf -t exon -g gene_id -o $filepath/countedfiles/$out_name $file
 
done
