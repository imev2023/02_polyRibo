#!/bin/bash

#BSUB -M 300M # Memory
#BSUB -R "select[mem>300M] rusage[mem=300M] span[hosts=1]"
#BSUB -n 1 # CPU's - matches hosts
#BSUB -q long # Queue (n=normal, l=long)
#BSUB -J nv4htseq
#BSUB -o /lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/023_Results/Logs/02_polyRibo_htseq_qlong_syes.log
#BSUB -e /lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/023_Results/Logs/02_polyRibo_htseq_qlong_syes.err

# run this from nv4/ --> bsub -G team267-grp < ./02_polyRibo/022_Scripts/02_polyRibo_htseq.sh

# load modules
module load htseq/2.0.3--py310h5aa3a86_1

# set environmental variables
WDIR="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/"
DATADIR_POLY="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/021_Data/fromKevin/data/polya/"
DATADIR_RIBO="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/021_Data/fromKevin/data/riboTest/"
RESDIR_RIBO="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/023_Results/02_polyRibo_htseq_riboOutput/"
RESDIR_POLY="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/023_Results/02_polyRibo_htseq_polyOutput/"
GTF="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/021_Data/Canis_lupus_familiaris.CanFam3.1.104.gtf"

for riboFile in $DATADIR_RIBO/*.cram
do
    riboName=$(basename $riboFile .cram) # Extract file name, without .cram. 
    echo "Processing $riboName"
    htseq-count -s yes -c "$RESDIR_RIBO/02_polyRibo_htseq_syes_$riboName.tsv" --additional-attr=gene_name $riboFile $GTF
done

# for polyFile in $DATADIR_RIBO/*.cram
# do
#     polyName=$(basename $polyFile .cram) # Extract file name, without .cram. 
#     echo "Processing $polyName"
#     htseq-count -s reverse -c "$RESDIR_RIBO/02_polyRibo_htseq_$polyName.tsv" --additional-attr=gene_name $polyFile $GTF
# done