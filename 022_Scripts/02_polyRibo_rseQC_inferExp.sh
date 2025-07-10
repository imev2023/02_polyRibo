#!/bin/bash

#BSUB -M 3.5GB # Memory
#BSUB -R "select[mem>3.5GB] rusage[mem=3.5GB] span[hosts=1]"
#BSUB -n 4 # CPU's
#BSUB -q normal # Queue (n=normal, l=long)
#BSUB -J nv4infExp
#BSUB -o /lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/023_Results/Logs/02_polyRibo_rseQC_inferExp_qNormal_m3.5GB_n4.log
#BSUB -e /lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/023_Results/Logs/02_polyRibo_rseQC_inferExp_qNormal_m3.5GB_n4.err

# run this from nv4/ --> bsub -G team267-grp < ./02_polyRibo/022_Scripts/02_polyRibo_rseQC_inferExp.sh

# load modules
module load subread/2.0.6--he4a0461_2
module load samtools/1.20--h50ea8bc_0
module load rseqc/4.0.0--py37h8902056_2
module load bedops/2.4.41--h4ac6f70_2

# set environmental variables. taking 2 samples each from Kevin's data to get an idea of 1) strandedness 2) using RNAseqQC for checking alignment stats. 
# 1) following https://bioinformaticsworkbook.org/dataAnalysis/RNA-Seq/Identify_Strandedness_Information.html#gsc.tab=0
# 2) following https://www.nature.com/articles/s41597-020-00719-4 (under methods, aligned reads distribution). Use geneBody_coverage for gene body cov of mapped reads. read_distribution for mapped reads distribute over a feature. tin.py is for RNA integrity at transcript level. 

WDIR="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/"
DATADIR_POLY="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/021_Data/fromKevin/test/polya/"
DATADIR_RIBO="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/021_Data/fromKevin/test/ribo/"
RESDIR_RIBO="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/023_Results/test/02_polyRibo_rseQC_inferExp_riboOutput/"
RESDIR_POLY="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/023_Results/test/02_polyRibo_rseQC_inferExp_polyOutput/"

# # Confirmed with Kevin. 
# GTF="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/021_Data/Canis_lupus_familiaris.CanFam3.1.104.gtf"

# # Convert gtf to bed file using bedops. Required for infer_experiment.py

# gtf2bed < $GTF > "/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/021_Data/Canis_lupus_familiaris.CanFam3.1.104.gtf.bed"

BED="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/021_Data/Canis_lupus_familiaris.CanFam3.1.104.gtf.bed"

for riboFile in $DATADIR_RIBO/*.cram
do
    # Extract file name, without .cram. 
    riboName=$(basename $riboFile .cram) 

    # Set temporary subdirectory per cram file. 
    TEMP_RESDIR_RIBO="$RESDIR_RIBO/$riboName/"

    # Create a subdirectory per cram file, this will not overwrite any existing directories.
    mkdir $TEMP_RESDIR_RIBO
    echo "Processing $riboName"
    # Bam files take up more space so creating and deleting within the loop.  
    samtools sort -@ 4 -m 500M -O bam -o "$DATADIR_RIBO/$riboName.bam" $riboFile 

    # Infer strandedness 
    infer_experiment.py --i "$DATADIR_RIBO/$riboName.bam" -r $BED

    # geneBody_coverage.py 

    # read_distribution.py

    # tin.py 

    # remove .bam file intermediate for space. 
    # rm -f "$DATADIR_RIBO/$riboName.bam"

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
    # Bam files take up more space so creating and deleting within the loop.  
    samtools sort -@ 4 -m 500M -O bam -o "$DATADIR_POLY/$polyName.bam" $polyFile 

    # Infer strandedness. Prints to err file
    infer_experiment.py --i "$DATADIR_POLY/$polyName.bam" -r $BED

    # geneBody_coverage.py 

    # read_distribution.py

    # tin.py 

    # remove .bam file intermediate for space. 
    # rm -f "$DATADIR_POLY/$polyName.bam"

done
