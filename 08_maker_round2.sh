# some of this code modified from https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2

# prepare gff files for round 2 of maker

cd /lustre/scratch/jmanthey/00_wbnu_maker/wbnu.maker.output

gff3_merge -n -s -d wbnu_master_datastore_index.log > wbnu_round1.all.maker.noseq.gff

# protein alignments
awk '{ if ($2 == "protein2genome") print $0 }' wbnu_round1.all.maker.noseq.gff > wbnu_round1.all.maker.protein2genome.gff
# repeat alignments
awk '{ if ($2 ~ "repeat") print $0 }' wbnu_round1.all.maker.noseq.gff > wbnu_round1.all.maker.repeats.gff

# move back up to parent directory
cd ..

# copy the original control file with options to save that information and modify the original control file
# for round 2
cp maker_opts.ctl maker_round1_opts.ctl

# rename the maker round 1 output directory
mv wbnu.maker.output/ maker_output_round1/

# modify the control file (copied in this directory)

# run maker round 2 (this took ~ 12 hours)
# script:

#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=maker
#SBATCH --partition quanah
#SBATCH --nodes=3 --ntasks=108
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=5G

source activate busco

mpiexec -n 106 maker



# summarize maker output

cd /lustre/scratch/jmanthey/00_wbnu_maker/wbnu.maker.output

gff3_merge -s -d wbnu_master_datastore_index.log > wbnu_round2.all.maker.gff
fasta_merge -d wbnu_master_datastore_index.log
gff3_merge -s  -n -d wbnu_master_datastore_index.log > wbnu_round2.all.maker.noseqs.gff
