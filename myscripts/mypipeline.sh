#!/bin/bash

# We'll use this a wrapper to run our different mapping scripts

# Path to my repo
myrepo="/users/b/w/bwworthi/Ecological-Genomics"

# echo ${myrepo}  <- to call my repo

# My population:
mypop="XSK"

# Directory to our cleaned and paired reads:
input="/data/project_data/RS_ExomeSeq/fastq/edge_fastq/pairedcleanreads/${mypop}"

# Directory to store the outputs of our mapping
# This part below goes to the common space
output="/data/project_data/RS_ExomeSeq/mapping"

# Run mapping.sh
source ./mapping.sh

# Run the post-processing steps
source ./process_bam.sh