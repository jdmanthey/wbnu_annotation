# move to working directory
cd /lustre/scratch/jmanthey/00_wbnu_maker

# make the control files 
maker -CTL

# edit the control files as necessary
# here, I used three protein files for the first round of maker:
# GCF_001522545.3_Parus_major1.1_protein.faa
# GCF_000247815.1_FicAlb1.5_protein.faa
# GCF_003957565.2_bTaeGut1.4.pri_protein.faa

# run maker round 1
# job script:

#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=maker
#SBATCH --partition quanah
#SBATCH --nodes=3 --ntasks=108
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=5G

mpiexec -n 106 maker

# job took longer than 48 hours, needed to restart (10 Dec 2021)
