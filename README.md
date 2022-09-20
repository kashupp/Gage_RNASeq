# Gage_RNASeq analysis pipeline.
Gage_RNASeq  is a  bioinformatics pipeline that is used to analyse RNA sequencing data obtained from genetically modified human cell lines. However, it can be used as general pupose RNASeq analysis pipeline. 

This pipeline is written using shell scripting and it can be run on UPPMAX using a single file, see run_RNASeq.sh.

This pipeline use FASTQC, MultiQC, Trimmomatic, STAR aligner and featureCounts for counting reads in genes.

Finally, differential gene expression analysis is performed using DSEq2 in R.
