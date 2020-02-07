import os
import re
import sys
import tempfile

snake_dir = os.path.dirname(workflow.snakefile) # workflow is a class in snakemake that point you to

# Taking care of environments
shell.executable("/bin/bash")
shell.prefix("source {}/env.sh; set -eo pipefail;".format(snake_dir)) # this line shuts down the program if any of the step fails
# -e shut down any script failed, -o shut down any pipe that failed 

# Load config file
if(os.path.exists("config.json")):
    cfile = "config.json"
else:
    cfile = snake_dir + "/config.json"
configfile:
    cfile


#
# Stolen from Michtell
# A little complicated to find the temp dir
# temp dir is a place programs can put temp files. 
# You will never have to change or think about anything here
#
''''
SSD_TMP_DIR = "/data/scratch/ssd"
if "TMPDIR" in os.environ:
    TMPDIR = os.environ['TMPDIR']
elif "TMPDIR" in config:
    TMPDIR = config['TMPDIR']
elif os.path.exists(SSD_TMP_DIR):
    TMPDIR = SSD_TMP_DIR
else:
    TMPDIR = tempfile.gettempdir()
'''


# Define variables
addSMtoheader = snake_dir + "/addSMtoheader.sh"
REGION_FILE = config["REGIONS"]
VCF = config["VCF"]
FOFN = config["FOFN"]
SAMPLE = config["SAMPLE"]
MINCOVERAGE = config["MINCOVERAGE"]
REF = config["REF"]
METHOD = config["METHOD"].lower()
SAM_FLAG = config["FLAG"]
PAD = config["PADDING"]
GSIZE = config["GENOME_SIZE"]

#
# Making sure all the regions are in samtools readable format
# 
def get_region():
	REGIONS = []
	with open(REGION_FILE, "r") as handle:
		for line in handle:
			line = re.split("[ \:\-\s\n_]", line)
			line = "{}_{}_{}".format(line[0], line[1], line[2])
			REGIONS.append(line)
	return REGIONS

wildcard_constraints:
	# whatshap only run on diploid so sorry banana
    # apparently some fern has 1260 chromosomes and hermit crabs
    # ,out of all the creatures, have 254 chromosomes
    # email me when you are running this on a hermit crab genome 
    region = "chr[0-9XYxy]+_\d+_\d+",
    haplotype = "h[12]",



#
#
# SNAKEMAKE
#
#

localrules:all

rule all:
	input:
		expand("asm/{region}/{haplotype}/pipeline.done", haplotype = ["h1","h2"], region = get_region())

#
# In order to run samtools merge or view with random accessing you
# need to index the bam file
# The same is needed for vcf files
#

rule index_input:
	input:
		vcf = VCF,
	output: VCF + ".tbi",
	resources:
		mem = 2,
	threads:1,
	params:
		log = "",
	shell:
		""" echo "$JOB_ID"
		tabix {input.vcf}
		"""

#
# extract reads and SNVs in a padded region
# all the input bams should be sorted and generated using PacBio aligners
#

rule extract_reads_and_SNVs:
	input:
		fofn = FOFN,
		vcf = VCF,
		tbi = VCF + ".tbi",
	output: 
		bam = temp("asm/{region}/compbined_reads.bam"),
		vcf = "asm/{region}/phased.vcf",
	resources:
		mem = 4,
	threads:1,
	params:
		log = "{region}",
	run:
		shell("echo $JOB_ID")
		region_lis = re.split("[\:\-\n_]", str(wildcards.region))
		region_lis[1] = max(0, int(region_lis[1]) - PAD)
		region_lis[2] = int(region_lis[2]) + PAD
		padded_region = "{}:{}-{}".format(region_lis[0], region_lis[1], region_lis[2])
		print("search region:" + padded_region)
		shell("mkdir asm/{} || true".format(wildcards.region))
		shell("samtools merge -R {} {} $(cat {})".format(padded_region, output.bam, input.fofn))
		shell("tabix {} {} > {}".format(input.vcf, padded_region, output.vcf))


#
# Whatshap requires the input files to be indexed
# so we have to do it for each region
#

rule index_region_files:
	input:
		bam = "asm/{region}/compbined_reads.bam",
		vcf = "asm/{region}/phased.vcf",
	output:
		bai = temp("asm/{region}/compbined_reads.bam.bai"),
		vcf = "asm/{region}/phased.vcf.gz",
		tbi = "asm/{region}/phased.vcf.gz.tbi",
	resources:
		mem = 4,
	threads:1,
	params:
		log = "{region}",
	shell:
		"""
		echo "$JOB_ID"
		cd asm/{wildcards.region}/
		# handle bam first
		samtools index compbined_reads.bam
		# then vcf 
		bgzip phased.vcf
		tabix phased.vcf.gz
		cd ../..
		"""


# As of whatshap 0.18, --ignore-read-group is not working on
# haplotag command. I have too add SM tag to the bam header
# this script add SM tag to @RG header
rule fix_header:
	input:
		bam = "asm/{region}/compbined_reads.bam",
	output:
		bai = temp("asm/{region}/compbined_reads.modified.bam.bai"),
		bam = temp("asm/{region}/compbined_reads.modified.bam"),
	resources:
		mem = 4,
	threads:1
	params:
		log = "{region}",
	shell:
		"""
		echo "$JOB_ID"
		{addSMtoheader} {SAMPLE} {input.bam} > {output.bam}
		samtools index {output.bam}
		"""


# Tag each read in the input bam file with HP:i:1 or HP:i:2
# based on the information provided by the phased VCF file
rule tag_read_by_haplotype:
	input:
		vcf = "asm/{region}/phased.vcf.gz",
		tbi = "asm/{region}/phased.vcf.gz.tbi",
		bai = "asm/{region}/compbined_reads.modified.bam.bai",
		bam = "asm/{region}/compbined_reads.modified.bam",
	output:
		bam = "asm/{region}/tagged.bam",
		bai = "asm/{region}/tagged.bam.bai",
	resources:
		mem = 4,
	threads: 1
	params:
		log = "{region}",
	shell:
		"""
		echo "$JOB_ID"
		if [[ $(stat {input.vcf} --printf "%s" | tr -d "\n") -lt 50 ]]; then
			echo "Pipeline Exception" >&2
			echo "Empty vcf file" >&2
			touch {output.bam}
			touch {output.bai}
		else
			whatshap haplotag \
			--ignore-read-groups \
			-o {output.bam} \
			--reference {REF} \
			{input.vcf} {input.bam}
			samtools index {output.bam}
		fi
		"""


# based on the haplotagging, split the tagged bam file into
# 2 bam files
rule partition_reads:
	input: "asm/{region}/tagged.bam",
	output:
		tmp = temp("asm/{region}/{haplotype}.tmp.sam"),
		bam = "asm/{region}/{haplotype}.bam",
	resources:
		mem = 4,
	threads:1
	params:
		log = "{region}",
	shell:
		"""
		echo "$JOB_ID"
		if ! [ -s {input} ]; then
			echo "Pipeline Exception" >&2
			echo "Empty tagged bam file" >&2
			touch {output.bam}
			touch {output.tmp}
		else
			samtools view -H {input} > {output.tmp}
			samtools view {input} | grep "HP:i:$(echo '{wildcards.haplotype}' | cut -c2-2)" >> {output.tmp} || true
			samtools view -bh {output.tmp} > {output.bam}
			echo "Read partitioning is complete" 
			count=$(samtools view {output.bam} -c)
			echo "$count reads are mapped to this region for {wildcards.haplotype}:"
			bedtools bamtobed -i {output.bam} | awk '{{print $1,$2,$3}}'
		fi
		"""

#
# Canu takes fasta or fastaq file. I am extracting sequences from
# bam files for each of the haplotype
# to override a error use command || true
# In bash == and = are the same and they perform string comparison
#

rule extract_sequence:
	input: "asm/{region}/{haplotype}.bam",
	output: temp("asm/{region}/{haplotype}.fasta"),
	resources:
		mem = 4,
	threads:1
	params:
		log = "{region}",
	shell:
		"""
		echo "$JOB_ID"
		if ! [ -s {input} ]; then
			echo "Pipeline Exception" >&2
			echo "Empty partitioned BAM file" >&2
			touch {output} 
		else
			samtools fasta {input} > {output} || true
		fi
		"""




#
# Assemble each region, 
# If the fasta input has less 10 reads this rule is skipped
# OMG this thing took me 5 hours to figure out, grep return 0 if there is no
# instance found !!!!
#

#
# all the condition checks that happen outside of a rule are run at the beginning of the pipeline
# which means that only user input can be used for conditional checks, you must handle empty output
# within a rule
#


rule assemble:
	input: 
		fasta = "asm/{region}/{haplotype}.fasta"
	output: 
		fasta = "asm/{region}/{haplotype}/asm.contigs.fasta",
		fai = "asm/{region}/{haplotype}/asm.contigs.fasta.fai",
		tiginfo = "asm/{region}/{haplotype}/asm.contigs.layout.tigInfo"
	resources:
		mem = 15,
	threads:4
	params:
		log = "{region}",
	run:
		shell('echo $JOB_ID')
		outfasta = output['fasta']
		outfai = output['fai']
		outtig = output["tiginfo"]
		region = wildcards.region
		haplotype = wildcards.haplotype
		infasta = input['fasta']

		if os.path.getsize(input['fasta']) == 0:
			shell(f"""
					echo "Pipeline Exception" >&2
					echo "Empty fasta file" >&2
					touch {outfasta}
					touch {outfai}
					touch {outtig}
				""")
		else:



			if METHOD in ['hifi']:
				canu_p = 'corMinCoverage=0 contigFilter="1 500 1.0 .25 1" \
					ovlMerThreshold=75 batOptions="-eg 0.01 -eM 0.01 -dg 6 -db \
					6 -dr 1 -ca 50 -cp 5"  correctedErrorRate=0.015 -pacbio-corrected'
				

				shell(f'''
			

					mkdir asm/{region}/{haplotype} || true
					mkdir $TMPDIR/asm || true
					cp {infasta} $TMPDIR
					pushd $TMPDIR
			

					canu {canu_p} {wildcards.haplotype}.fasta \
						stopOnLowCoverage=1 \
						genomeSize={GSIZE} \
						corOutCoverage=300 \
						corMhapSensitivity=high \
						-p asm useGrid=false  \
						-d asm \
						maxThreads=4 || true
			


					samtools faidx asm/asm.contigs.fasta || \
						(touch asm/asm.contigs.fasta && \
						touch asm/asm.contigs.fasta.fai)
					popd
					mv $TMPDIR/asm/asm.contigs.fasta asm/{region}/{haplotype} || touch asm/{region}/{haplotype}/asm.contigs.fasta
					mv $TMPDIR/asm/asm.contigs.fasta.fai asm/{region}/{haplotype} || touch asm/{region}/{haplotype}/asm.contigs.fasta.fai
					mv $TMPDIR/asm/asm.contigs.layout.tigInfo asm/{region}/{haplotype} || touch asm/{region}/{haplotype}/asm.contigs.layout.tigInfo
					cp $TMPDIR/asm/canu-logs/*_canu asm/{region}/{haplotype}
					''')

			




			elif METHOD in ['clr']:
				canu_p	= 'corMinCoverage=1 contigFilter="1 500 1.0 .75 1" \
				cnsThreads=1 ovlThreads=1 gnuplot=undef'
				shell(f"""
					mkdir asm/{region}/{haplotype} || true
					mkdir $TMPDIR/asm || true
					cp {infasta} $TMPDIR
					pushd $TMPDIR
					canu {canu_p} -pacbio-raw {wildcards.haplotype}.fasta \
stopOnLowCoverage=1 \
genomeSize={GSIZE} \
corMhapSensitivity=high \
-p asm useGrid=false  \
-d asm \
mhapThreads=4 || true	
					samtools faidx asm/asm.contigs.fasta || \
							(touch asm/asm.contigs.fasta && \
							touch asm/asm.contigs.fasta.fai)
					popd
					mv $TMPDIR/asm/asm.contigs.fasta asm/{region}/{haplotype} || touch asm/{region}/{haplotype}/asm.contigs.fasta
					mv $TMPDIR/asm/asm.contigs.fasta.fai asm/{region}/{haplotype} || touch asm/{region}/{haplotype}/asm.contigs.fasta.fai
					mv $TMPDIR/asm/asm.contigs.layout.tigInfo asm/{region}/{haplotype} || touch asm/{region}/{haplotype}/asm.contigs.layout.tigInfo
					cp $TMPDIR/asm/canu-logs/*_canu asm/{region}/{haplotype}
					""")


				
			
			else:
				shell('echo "Method not supported. Please modify the method variable in the config file."')
				shell('echo "Supported methods: Hifi, CLR."')
				sys.exit(1)

#
# You have to map the contigs back to the referance
#
#


rule map_contigs:
	input:
		bam = "asm/{region}/{haplotype}.bam",
		fasta = "asm/{region}/{haplotype}/asm.contigs.fasta",
		fai = "asm/{region}/{haplotype}/asm.contigs.fasta.fai",
	output:
		bam = temp("asm/{region}/{haplotype}/asm.blazr.bam"),
		pbi = temp("asm/{region}/{haplotype}/asm.blazr.bam.pbi"),
	resources:
		mem = 4,
	threads:1
	params:
		log = "{region}",
	shell:
		"""
		echo "$JOB_ID"
		if ! [ -s {input.fasta} ]; then
			echo "Pipeline Exception" >&2
			echo "Empty assembled sequence" >&2
			touch {output.bam}
			touch {output.pbi}
		elif ! [ -s {input.bam} ]; then
			echo "Pipeline Exception" >&2
			echo "Empty partitioned BAM file" >&2
			touch {output.bam}
			touch {output.pbi}
		else
			pbmm2 align --preset SUBREAD -j 4 {input.fasta} {input.bam} | samtools view -u -F {SAM_FLAG} - | samtools sort - > {output.bam}
			samtools index {output.bam} && pbindex {output.bam}
		fi
		"""

#
# Do Arrow polishing on the contigs
#

rule polish:
	input: 
		bam = "asm/{region}/{haplotype}/asm.blazr.bam",
		pbi = "asm/{region}/{haplotype}/asm.blazr.bam.pbi",
		fai = "asm/{region}/{haplotype}/asm.contigs.fasta.fai",
		fasta = 'asm/{region}/{haplotype}/asm.contigs.fasta',

	output: 
		fasta = "asm/{region}/{haplotype}/asm.corrected.fasta",
		sam = temp("asm/{region}/{haplotype}/temp.contig.sam"),
		fastq = temp("asm/{region}/{haplotype}/temp.reads.fastq"),

	resources:
		mem = 4,
	threads:1
	params:
		log = "{region}",
	run:
		shell('echo $JOB_ID')
		
		infasta = input['fasta']
		outfasta = output['fasta']
		fastq = output['fastq']
		bam = input['bam']
		sam = output['sam']

		if os.path.getsize(input['fasta']) == 0:
			shell(f'''
				echo "Pipeline Exception" >&2
				echo "Empty blazr mapping" >&2
				touch {outfasta} {sam} {fastq}			
				''')

		elif METHOD in ['hifi']:
			shell(f'''
					echo "Racon polishing"
					samtools view {bam} > {sam}"
					samtools fastq {bam} > {fastq}"
				''')

			shell(f'''
					racon {fastq} {sam} {infasta} -u -t {threads} > {outfasta}
				''')

		elif METHOD in ['clr']:
			shell(f'''
					echo "Arrow polishing"
					touch {sam}
					touch {fastq}
				''')
			shell(f'''
					arrow --noEvidenceConsensusCall nocall \
					--minCoverage 3 -j 1 -r {infasta} \
					-o {outfasta} {bam}
				''')


rule stats:
	input: "asm/{region}/{haplotype}/asm.corrected.fasta"
	output: "asm/{region}/{haplotype}/pipeline.done"
	resources:
		mem = 1,
	threads:1
	params:
		log = "{region}"
	shell:
		'''
		echo "===============================" > {output}
		echo "Region {wildcards.region}" >> {output}
		echo {METHOD} assembly
		homo_count=$(less asm/{wildcards.region}/phased.vcf.gz | grep -o -P "1\|1" | wc -l) || true
		h1_count=$(less asm/{wildcards.region}/phased.vcf.gz | grep -o -P "0\|1" | wc -l) || true
		h2_count=$(less asm/{wildcards.region}/phased.vcf.gz | grep -o -P "1\|0" | wc -l) || true
		read_count=$(samtools view -c asm/{wildcards.region}/compbined_reads.bam) || true
		echo "unphased SNVs: $homo_count" >> {output}
		echo "phased SNVs: $h1_count + $h2_count" >> {output}
		echo "h1: $h1_count" >> {output}
		echo "h2: $h2_count" >> {output}
		echo "reads: $read_count" >> {output}
		echo "===============================" >> {output}
		partitioned_read_count=$(samtools view -c asm/{wildcards.region}/{wildcards.haplotype}.bam) || true
		echo "Haplotype {wildcards.haplotype}" >> {output}
		if [ -s "asm/{wildcards.region}/{wildcards.haplotype}/asm.corrected.fasta" ]; then
			echo "Assembly done" >> {output}
			echo "Coverage:" >> {output}
			awk '(NR!=1 && $4 != 0.00){{print $5}}' asm/{wildcards.region}/{wildcards.haplotype}/asm.contigs.layout.tigInfo >> {output} || true
		else
			echo "Assembly failed" >> {output}
		fi
		echo "Partitioned reads: $partitioned_read_count" >> {output}
		echo "===============================" >> {output}
		echo "Debug"
		echo "{input}"
		'''










"""

The following rules are deprecated
Checkpoint functionality is very difficult to use when handling potential error by the pipeline


rule termination:
	input: "asm/{region}/{haplotype}.fasta"
	output: "asm/{region}/{haplotype}.termination_signal"
	resources:
		mem = 1,
	threads:1
	params:
		log = "{region}.{haplotype}",
	shell:
		'''
		echo "$JOB_ID"
		touch {output}
		'''

def check_assembly(wildcards):
	fasta = checkpoints.assemble.get(**wildcards).output[0]
	if os.path.exists(fasta) and os.path.getsize(fasta):
		return "asm/{region}/{haplotype}/asm.arrow_polished.fasta"
	else:
		return "asm/{region}/{haplotype}.termination_signal"


			blasr {input.bam} {input.fasta} \
--clipping subread --bam --bestn 1 \
--nproc 1 --out - | samtools sort -m 4G -T tmp -o {output.bam}
			pbindex {output.bam}
"""
'''

		echo "$JOB_ID"
		if ! [ -s {input.bam} ]; then
			echo "Pipeline Exception" >&2
			echo "Empty blazr mapping" >&2
			touch {output}
		elif [[ $(echo {METHOD} | tr '[:upper:]' '[:lower:]') == 'hifi' ]]; then
			echo "RACON polishing"
			racon {output.tmp_fastq} {output.tmp_sam} {input.asm} -u -t {threads} > {output.tmp_cor} || \
				( >&2 echo " RACON FAILED TO RUN " && cp {input.asm} {output.tmp_cor} ) """)	
		else
		fi
		"""
'''