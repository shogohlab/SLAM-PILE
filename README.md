# SLAM-PILE v0.1
SLAM-PILE is a T > C alignment software package based on STAR and Rsamtools to recover T > C conversions from SLAMseq datasets. SLAM-PILE analyzes T > C alignments covering the entire mRNA length instead of just the 3’ end of mRNA.

SLAM-PILE requires the following dependencies: 
Rsamtools, pandas, samtools.

SLAM-PILE also requires:
A) input bam file generated from mapping catabolic SLAMseq fastq data using 
STAR 2.7.10a --outMultimapperOrder Random --outSAMmultNmax 1
B) bam index file for the above bam file generated using
samtools index


1)
Denote a  "destination" folder, gather requisite bam and bai files in this folder then run
SLAM-PILE1.R

ouput: ".pileup.csv" file

2) 
Edit the following within SLAM-PILE2.py
a) gtf file
b) genome.fa file
c)  filename of requisite ".pileup.csv" file

then run
SLAM-PILE2.py

ouput:  ".conversionStatistic.csv" file.

3)
Use preferred curve-fitting process to calculate gene half-life using ".conversionStatistic.csv" file
