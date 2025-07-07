#!/bin/bash

#BSUB -M 50 # Memory
#BSUB -R "select[mem>50] rusage[mem=50] span[hosts=1]"
#BSUB -n 1 # CPU's - matches hosts
#BSUB -q normal # Queue (n=normal)
#BSUB -J nv4flagstats
#BSUB -o /lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/023_Results/Logs/02_polyRibo_flagstats.log
#BSUB -e /lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/023_Results/Logs/02_polyRibo_flagstats.err

# run this from the output dir bsub -G team267-grp < ../../022_Scripts/02_polyRibo_flagstats.sh

# load modules
module load samtools/1.20--h50ea8bc_0

# set environmental variables
WDIR="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/"
DATADIR_POLY="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/021_Data/fromKevin/data/polya/"
DATADIR_RIBO="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/021_Data/fromKevin/data/ribo/"
RESDIR_RIBO="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/023_Results/02_polyRibo_flagstats_riboOutput/"
RESDIR_POLY="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/023_Results/02_polyRibo_flagstats_polyOutput/"


for riboFile in $DATADIR_RIBO/*.cram
do
    riboName=$(basename $riboFile .cram) # Extract file name, without .cram. 
    echo "Processing $riboName"
    samtools flagstat $riboFile -O tsv > "$RESDIR_RIBO/02_polyRibo_flagstats_$riboName.tsv"
    # samtools stat $riboFile > "$RESDIR_RIBO/02_polyRibo_stats_$riboName.tsv"
    # samtools sort -O CRAM $riboFile > "$RESDIR_RIBO/02_polyRibo_sorted_$riboName.cram"
done

for polyFile in $DATADIR_POLY/*.cram
do
    polyName=$(basename $polyFile .cram) # Extract file name, without .cram. 
    echo "Processing $polyName"
    samtools flagstat $polyFile -O tsv > "$RESDIR_POLY/02_polyRibo_flagstats_$polyName.tsv"
    # samtools stat $polyFile > "$RESDIR_POLY/02_polyRibo_stats_$polyName.tsv"
done

# test sort on 1 file for htseq count downstream
# samtools sort -O CRAM -o "/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/023_Results/02_polyRibo_htseq_riboOutput/testSamtoolsSort.cram" "nv4/02_polyRibo/021_Data/fromKevin/data/ribo/2169Ta.ribo.cram"