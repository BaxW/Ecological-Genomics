#!/bin/bash

cd ~/Ecological-Genomics/myresults

mkdir fastqc

for file in /data/project_data/RS_ExomeSeq/fastq/edge_fastq/XSK*fastq.gz

do

fastqc ${file} -o fastqc/

done
