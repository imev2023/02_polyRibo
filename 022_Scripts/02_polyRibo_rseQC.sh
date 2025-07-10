#!/bin/bash
RESDIR="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/023_Results/02_polyRibo_rseQC_riboOutput/"

# # Confirmed with Kevin. 
# GTF="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/021_Data/Canis_lupus_familiaris.CanFam3.1.104.gtf"

# # Convert gtf to bed file using bedops. Required for infer_experiment.py. Not a BED12, not sure if this is a problem. 
BED_GAPDH="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/020_Software/result.bed"

# This will be read in from samples.csv (see execute script)
file="$1" 

# Extract file name, without .bam. 
name=$(basename $file .bam) 

# Set temporary subdirectory per bam file. 
TEMP_RESDIR="$RESDIR/$name/"
mkdir $TEMP_RESDIR
echo "Processing $name"

# My bed file
## 1. 5' to 3' bias
geneBody_coverage.py -r $BED_GAPDH -i $file -o "${TEMP_RESDIR}_${name}_geneBody_BED_GAPDH"

## 2. Intron vs exon origin
infer_experiment.py -r $BED_GAPDH -i $file > "${TEMP_RESDIR}_${name}_infer_BED_GAPDH.txt"

## 3. Read distribution by feature type
read_distribution.py -r $BED_GAPDH -i $file > "${TEMP_RESDIR}_${name}_readDist_BED_GAPDH.txt"