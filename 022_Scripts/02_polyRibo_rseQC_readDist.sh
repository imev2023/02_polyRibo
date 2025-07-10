#!/bin/bash

#BSUB -M 3.5GB # Memory
#BSUB -R "select[mem>3.5GB] rusage[mem=3.5GB] span[hosts=1]"
#BSUB -n 1 # CPU's
#BSUB -q normal # Queue (n=normal, l=long)
#BSUB -J nv4readDist
#BSUB -o /lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/023_Results/Logs/02_polyRibo_rseQC_readDist_qNormal_m3.5GB_n1.log
#BSUB -e /lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/023_Results/Logs/02_polyRibo_rseQC_readDist_qNormal_m3.5GB_n1.err

# run this from nv4/ --> bsub -G team267-grp < ./02_polyRibo/022_Scripts/02_polyRibo_rseQC_readDist.sh

# load modules
# module load subread/2.0.6--he4a0461_2
module load samtools/1.20--h50ea8bc_0
module load rseqc/4.0.0--py37h8902056_2
# module load bedops/2.4.41--h4ac6f70_2

# set environmental variables. taking 2 samples each from Kevin's data to get an idea of 1) strandedness 2) using RNAseqQC for checking alignment stats. 
# 1) following https://bioinformaticsworkbook.org/dataAnalysis/RNA-Seq/Identify_Strandedness_Information.html#gsc.tab=0
# 2) following https://www.nature.com/articles/s41597-020-00719-4 (under methods, aligned reads distribution). Use geneBody_coverage for gene body cov of mapped reads. read_distribution for mapped reads distribute over a feature. tin.py is for RNA integrity at transcript level. 

WDIR="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/"
DATADIR_POLY="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/021_Data/fromKevin/test/polya/"
DATADIR_RIBO="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/021_Data/fromKevin/test/ribo/"
RESDIR_RIBO="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/023_Results/test/02_polyRibo_rseQC_readDist_riboOutput/"
RESDIR_POLY="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/023_Results/test/02_polyRibo_rseQC_readDist_polyOutput/"

# # Confirmed with Kevin. 
# GTF="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/021_Data/Canis_lupus_familiaris.CanFam3.1.104.gtf"

# # Convert gtf to bed file using bedops. Required for infer_experiment.py. Not a BED12, not sure if this is a problem. 

# gtf2bed < $GTF > "/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/021_Data/Canis_lupus_familiaris.CanFam3.1.104.gtf.bed"

BED="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/021_Data/Canis_lupus_familiaris.CanFam3.1.104.gtf.bed12"

# Index the bam files. 
# samtools index -@ 2 -M "/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/021_Data/fromKevin/test/rseQC/2169Ta.polya.bam" "/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/021_Data/fromKevin/test/rseQC/2169Ta.ribo.bam"

# To compute gene body cov, want to check on all bed file, but also on subset of housekeeping genes (like in RSEQC tutorial --> although maybe we can filter on these?)
# For the test data, I had already computed the bam files, but would need to do this for all files if we're interested. 
# Use xargs?
# -i input directory with more than one bam file (must be sorted and indexed). For larger group will use option 4 (txt file with bamfile per row.)
read_distribution.py -i "/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/021_Data/fromKevin/test/rseQC/2169Ta.polya.bam" \
 -r $BED > "/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/023_Results/test/02_polyRibo_rseQC_readDist_test/02_polyRibo_rseQC_readDist_test_2169Ta.polya.bam_result.txt"

sed -n "5,15p" "/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/023_Results/test/02_polyRibo_rseQC_readDist_test/02_polyRibo_rseQC_readDist_test_2169Ta.polya.bam_result.txt" > "/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/023_Results/test/02_polyRibo_rseQC_readDist_test/02_polyRibo_rseQC_readDist_test_2169Ta.polya.bam_result.UTR.txt"