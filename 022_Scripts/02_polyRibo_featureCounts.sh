#!/bin/bash

#BSUB -M 3.5GB # Memory
#BSUB -R "select[mem>3.5GB] rusage[mem=3.5GB] span[hosts=1]"
#BSUB -n 4 # CPU's
#BSUB -q normal # Queue (n=normal, l=long)
#BSUB -J nv4fC
#BSUB -o /lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/023_Results/Logs/02_polyRibo_featureCounts_qLong_m3.5G_n4.log
#BSUB -e /lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/023_Results/Logs/02_polyRibo_featureCounts_qLong_m3.5G_n4.err

# run this from nv4/ --> bsub -G team267-grp < ./02_polyRibo/022_Scripts/02_polyRibo_featureCounts.sh

# load modules
module load subread/2.0.6--he4a0461_2
module load samtools/1.20--h50ea8bc_0

# set environmental variables
WDIR="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/"
DATADIR_POLY="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/021_Data/fromKevin/data/polya/"
DATADIR_RIBO="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/021_Data/fromKevin/data/ribo/"
RESDIR_RIBO="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/023_Results/02_polyRibo_featureCounts_riboOutput/"
RESDIR_POLY="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/023_Results/02_polyRibo_featureCounts_polyOutput/"

# Confirmed with Kevin. 
GTF="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/021_Data/Canis_lupus_familiaris.CanFam3.1.104.gtf"

for riboFile in $DATADIR_RIBO/*.cram
do
    # Extract file name, without .cram. 
    riboName=$(basename $riboFile .cram) 

    # Set temporary subdirectory per cram file. 
    TEMP_RESDIR_RIBO="$RESDIR_RIBO/$riboName/"

    # Create a subdirectory per cram file, this will not overwrite any existing directories.
    mkdir $TEMP_RESDIR_RIBO
    echo "Processing $riboName"
    # featureCounts is a lot faster than htseq, but requires bam rather than cram files. Bam files take up more space so creating and deleting within the loop.  
    samtools sort -@ 4 -m 500M -O bam -o "$DATADIR_RIBO/$riboName.bam" $riboFile #

    # run featureCounts two ways. This is assuming -stranded = reverse, following Maurine's seq. This may be incorrect. 
    featureCounts -T 4 -o "$TEMP_RESDIR_RIBO/02_polyRibo_featureCounts_s2_$riboName" -p -s 2 -a $GTF "$DATADIR_RIBO/$riboName.bam" # T = threads, -o = output, -s stranded (2 for reverse?), -p paired-end reads
    featureCounts -T 4 -f -o "$TEMP_RESDIR_RIBO/02_polyRibo_featureCounts_fs2_$riboName" -p -s 2 -a $GTF "$DATADIR_RIBO/$riboName.bam" # T = threads, -o = output, -s stranded (2 for reverse?), -f performs read counting at exon level rather than genes. 

    # remove .bam file intermediate for space. 
    rm -f "$DATADIR_RIBO/$riboName.bam"

done

for polyFile in $DATADIR_POLY/*.cram
do
    # Extract file name, without .cram. 
    polyName=$(basename $polyFile .cram) 

    # Set temporary subdirectory per cram file. 
    TEMP_RESDIR_POLY="$RESDIR_POLY/$polyName/"

    # Create a subdirectory per cram file, this will not overwrite any existing directories.
    mkdir $TEMP_RESDIR_POLY

    echo "Processing $polyName"
    # featureCounts is a lot faster than htseq, but requires bam rather than cram files. Bam files take up more space so creating and deleting within the loop.  
    samtools sort -@ 4 -m 500M -O bam -o "$DATADIR_POLY/$polyName.bam" $polyFile #

    # run featureCounts two ways. This is assuming -stranded = reverse, following Maurine's seq. Confirmed with Kevin. 
    # T = threads, -o = output, -s stranded (2 for reverse), -f performs read counting at exon level rather than genes. -R assignment results per read/read-pair, saved in -o. 
    featureCounts -T 4 -o "$TEMP_RESDIR_POLY/02_polyRibo_featureCounts_s2_$polyName" -p -s 2 --extraAttributes "gene_id" -a $GTF "$DATADIR_POLY/$polyName.bam"
    featureCounts -T 4 -f -o "$TEMP_RESDIR_POLY/02_polyRibo_featureCounts_fs2_$polyName" -p -s 2 --extraAttributes "gene_id" -a $GTF "$DATADIR_POLY/$polyName.bam" 

    # remove .bam file intermediate for space. 
    rm -f "$DATADIR_POLY/$polyName.bam"

done