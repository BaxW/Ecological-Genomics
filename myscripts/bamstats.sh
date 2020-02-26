#!/bin/bash

# set repo

myrepo="/users/b/w/bwworthi/Ecological-Genomics/"

mypop="XSK"

output="/data/project_data/RS_ExomeSeq/mapping/"

echo "Num.reads R1 R2 Paired MateMapped Singletons MateMappedDiffChr" > ${myrepo}/myresults/${mypop}.flagstats.txt

for file in ${output}/BWA/${mypop}*.sorted.rmdup.bam

	do 
	
f=${file/.sorted.rmdup.bam}
name=`basename ${f}`
echo $name >> ${myrepo}/myresults/${mypop}.names.txt

# take flagstat output, grab only rows between 6 and 12 and then take first col and flip it into a row 
samtools flagstat ${file} | awk 'NR>=6&&NR<=12 {print $1}' | column -x 

done >> ${myrepo}/myresults/${mypop}.flagstats.txt