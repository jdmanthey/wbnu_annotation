#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=RM
#SBATCH --partition quanah
#SBATCH --nodes=1 --ntasks=18
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --array=1-39

# use a custom repeatmasker database to annotate the genome for TEs 
# includes the:
# RepBase vertebrate database v24.03 sequences
# certhia americana custom repeatmodeler sequences (doi: 10.1093/gbe/evab120)
# colaptes auratus custom repeatmodeler sequences (doi:10.1093/g3journal/jkaa026)

# run repeat masker v1.332
cd /lustre/scratch/jmanthey/wbnu_scaffolds
RepeatMasker -pa 18 -s -lib ~/RepeatMasker/Libraries/custom_library_certhia_colaptes.fa Scaffold_${SLURM_ARRAY_TASK_ID}.fasta
