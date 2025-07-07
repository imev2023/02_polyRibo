#!/bin/bash

#BSUB -M 5GB # Memory
#BSUB -R "select[mem>5GB] rusage[mem=5GB] span[hosts=2]"
#BSUB -n 2 # CPU's - matches hosts
#BSUB -q normal # Queue (n=normal)
#BSUB -J nv4Flagstats
#BSUB -o /lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/023_Results/Logs/02_polyRibo_flagstats.log
#BSUB -e /lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/023_Results/Logs/02_polyRibo_flagstats.err

# run this from the output dir bsub -G team267-grp < ../../022_Scripts/02_polyRibo_flagstats.sh

# load modules
module load samtools/1.20--h50ea8bc_0

# set environmental variables
WDIR="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/"
POLYDIR="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/021_Data/fromKevin/data/polya/"
DATADIR_RIBO="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/021_Data/fromKevin/data/ribo/"
RESDIR_RIBO="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/023_Results/02_polyRibo_flagstats_riboOutput/"

for riboFile in $DATADIR_RIBO/*.cram
do
    riboName=$(basename $riboFile .cram) # Extract file name, without .cram. 
    echo "Processing $riboName"
    samtools flagstat $riboFile -O tsv > "$RESDIR_RIBO/02_polyRibo_flagstats_$riboName.tsv"
done