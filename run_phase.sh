# Getting the path where this script live
# Snake file should be in the directory as this file
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Taking care of the environment
if [ -s env.sh ]; then
    source env.sh
else
    source $DIR/env.sh
fi


# Define variables
SNAKEFILE=$DIR/phase.snakemake.py
JOBNUM=150
WAITTIME=100
LOGDIR=log
RETRY=0
mkdir -p $LOGDIR

snakemake -s $SNAKEFILE \
        --drmaa " -P eichlerlab \
                -q eichler-short.q \
                -l h_rt=48:00:00  \
                -l mfree={resources.mem}G \
                -l gpfsstate=0 \
                -pe serial {threads} \
                -V -cwd \
                -o $LOGDIR/{rule}.{params.log}.o \
                -e $LOGDIR/{rule}.{params.log}.e \
                -S /bin/bash" \
        --jobs $JOBNUM \
        --latency-wait $WAITTIME \
        --restart-times $RETRY \
		--rerun-incomplete \
	-p \
	-k $1 $2
