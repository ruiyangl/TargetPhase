#!/bin/bash
sample_name=$1
file=$2
# check if sample name exist in the file, if not  
if ! $(samtools view $file -H | grep -qi $sample_name -)
then
 	echo "No sample name found" >&2
	echo "The sample name given is" $sample_name >&2
	samtools view -h $file | awk -v sample_name=$sample_name 'BEGIN {OFS = "\t"} ($1!="@RG"){print} ($1=="@RG"){print $0, "SM:"sample_name}' | samtools view -bh -
else
	samtools view -b -h $file
fi
