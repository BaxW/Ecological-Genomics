#!/bin/bash

# This script will run the read mapping using "bwa" program
ref="/data/project_data/RS_ExomeSeq/ReferenceGenomes/Pabies1.0-genome_reduced.fa"

# Write a loop to map each individual within my population

for forward in ${input}*_R1.cl.pd.fq # * wild card to read everything that starts with my pop id

do
	reverse=${forward/_R1.cl.pd.fq/_R2.cl.pd.fq} # replace _R1.. with _R2..
	f=${forward/_R1.cl.pd.fq/} # replace _R1.. with nothing
	name=`basename ${f}`
	bwa mem -t 4 -M ${ref} ${forward} ${reverse} > ${output}/BWA/${name}.sam
done 
