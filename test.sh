js=$1

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
