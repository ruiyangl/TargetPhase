!#!/bin/bash

# Use this script to submit large region list
# split region file into multiple smaller batches and submit them sequentially
# remove split file after the batch is finished
# usage bash submit_large_catalog json_file do_split 





js=$1

find -maxdepth 1 -type f -name "splited_region*" | grep "."

if [[ $? -ne 0 ]]; then
region_file=$(
python << CHICKEN  ## get region file name
import json
import sys

with open("$js", 'r') as infile:
	data =  json.load(infile)
	region = data['REGIONS']
	sys.stdout.write(region)
CHICKEN
)


split -l 2000 $region_file splited_region
fi


find -maxdepth 1 -type f -name "splited_region*" | while read file; do export region=$file;
python << MOOOO
import json

with open("$js", 'r') as infile:
	data = json.load(infile)
data['REGIONS'] = "$region"
with open("$js", 'w') as outfile:
	json.dump(data, outfile)
MOOOO
/net/eichler/vol27/projects/ruiyang_projects/nobackups/vntr_project/phasedsv/Xphasing/run_phase_no_log.sh
if [[ $? -eq 0 ]]; then
	echo "Removing finished region file $region"
	rm $region
fi

# Cleaning up
echo "Snakedir clean up"
rm -r .snakemake/{conda-archive\
	,conda,singularity,shadow\
	tmp.d_lyyizx,metadata,locks}
rm -r .snakemake/tmp*

done
