#!/bin/bash

#BSUB -M 1000 # Memory
#BSUB -R "select[mem>1000] rusage[mem=1000] span[hosts=1]"
#BSUB -n 4 # CPU's
#BSUB -q normal # Queue (n=normal, l=long)
#BSUB -J nv4rseqcpoly
#BSUB -o /lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/023_Results/Logs/02_polyRibo_rseQC_polyOutput_execute_qNormal_m1000_n4.log
#BSUB -e /lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/023_Results/Logs/02_polyRibo_rseQC_polyOutput_execute_qNormal_m1000_n4.err

# run this from nv4/ --> bsub -G team267-grp < ./02_polyRibo/022_Scripts/02_polyRibo_rseQC_poly_execute.sh

# load modules
module load samtools/1.20--h50ea8bc_0
module load rseqc/4.0.0--py37h8902056_2
module load gffread/0.12.7--hdcf5f25_4 

# GTF_GAPDH="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/021_Data/02_polyARibo_DESEQ_ENSCAFG00000015077.gtf"
# GTF_FULL="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/021_Data/Canis_lupus_familiaris.CanFam3.1.104.gtf"
# # convert hskg GAPDH bed file --> just testing. 
# gffread $GTF_GAPDH -T -o- | gtfToGenePred - /dev/stdout | genePredToBed - 02_polyARibo_DESEQ_ENSCAFG00000015077_bed12_output.bed
# ./gtfToGenePred $GTF_GAPDH test.genePhred && ./genePredToBed test.genePhred result.bed && rm test.genePhred
# ./gtfToGenePred $GTF_FULL test.genePhred && ./genePredToBed test.genePhred Canis_lupus_familiaris.CanFam3.1.104_result.bed && rm test.genePhred

# make samples.csv
# Set this to your BAM folder (absolute or relative)
BAM_DIR_POLY="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/021_Data/bam/polya/"

# Output file
OUT_CSV="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/021_Data/bam/poly/02_polyRibo_rseQC_polyOutput_bam.csv"

# # Write header
echo "bam" > "$OUT_CSV"

# # Find all .bam files in the folder and append to CSV
find "$BAM_DIR_POLY" -type f -name "*.bam" | sort >> "$OUT_CSV"

# echo "Created $OUT_CSV with $(wc -l < $OUT_CSV) entries (including header)"

RSEQC_SCRIPT="/lustre/scratch125/casm/staging/team267_murchison/nv4/02_polyRibo/022_Scripts/02_polyRibo_rseQC_poly.sh"

chmod u+x $RSEQC_SCRIPT
# tail -n +2 samples.csv skips the header, -n1 reads one column per _samples.csv, -P4 runs with 4 CPUs. 
tail -n +2 $OUT_CSV | xargs -n1 -P4 $RSEQC_SCRIPT
