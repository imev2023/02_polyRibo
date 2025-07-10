#!/bin/bash

#BSUB -M 60M # Memory
#BSUB -R "select[mem>60M] rusage[mem=60MB] span[hosts=1]"
#BSUB -n 1 # CPU's
#BSUB -q normal # Queue (n=normal, l=long)
#BSUB -J nv4idxPoly
#BSUB -o /lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/023_Results/Logs/02_polyRibo_samtoolIdx_poly_qNormal_60M_n1.log
#BSUB -e /lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/023_Results/Logs/02_polyRibo_samtoolsIdx_poly_qNormal_60M_n1.log

# run this from nv4/ --> bsub -G team267-grp < ./02_polyRibo/022_Scripts/02_polyRibo_samtoolsSort_poly.sh

# load modules
module load subread/2.0.6--he4a0461_2
module load samtools/1.20--h50ea8bc_0

# set environmental variables
WDIR="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/"
DATADIR_POLY="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/021_Data/fromKevin/data/polya/"
DATADIR_RIBO="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/021_Data/fromKevin/data/ribo/"
RESDIR_RIBO="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/023_Results/02_polyRibo_samtoolsSort_riboOutput/"
RESDIR_POLY="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/023_Results/02_polyRibo_samtoolsSort_polyOutput/"
NFSDIR_RIBO="/nfs/users/nfs_n/nv4/02_polyRibo/023_Results/02_polyRibo_featureCounts_riboOutput/"
NFSDIR_POLY="/nfs/users/nfs_n/nv4/02_polyRibo/023_Results/02_polyRibo_featureCounts_polyOutput/"

NEW_DATADIR_POLY="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/021_Data/bam/polya/"

# for polyFile in $DATADIR_POLY/*.cram
# do
#     # Extract file name, without .cram. 
#     polyName=$(basename $polyFile .cram) 

#     # Set temporary subdirectory per cram file. 
#     TEMP_RESDIR_POLY="$RESDIR_POLY/$polyName/"
#     # Create a subdirectory per cram file, this will not overwrite any existing directories.
#     mkdir $TEMP_RESDIR_POLY


#     echo "Processing $polyName"
    
#     # featureCounts is a lot faster than htseq, but requires bam rather than cram files. Bam files take up more space so creating and deleting within the loop.  
#     samtools sort -@ 8 -m 2G -O bam -o "$NEW_DATADIR_POLY/$polyName.bam" $polyFile #

# done

for polyFile in $NEW_DATADIR_POLY/*.bam
do
   
   samtools index $polyFile

done
