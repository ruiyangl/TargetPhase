#!/bin/bash

module purge
/etc/profile.d/modules.sh
module load modules modules-init modules-gs/prod modules-eichler
module load gcc/8.1.0
module load bedtools/2.27.1 samtools/1.9\
 miniconda/4.5.12 whatshap/0.18\
 canu java/8u151  tabix/0.2.6 pbconda