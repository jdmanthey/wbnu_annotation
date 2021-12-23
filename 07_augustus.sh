# some parts of code from https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2

interactive -p nocona -c 128

cd /lustre/scratch/jmanthey/00_wbnu_maker/wbnu.maker.output

# get gff from maker without fasta seqs
gff3_merge -n -s -d wbnu_master_datastore_index.log > wbnu_rnd1.all.maker.noseq.gff

cd ..
mkdir -p augustus/round1
cd augustus/round1

# extract the fasta sequences and surrounding 1000 bp
# some will be missed because bedtools will error when the 1000 bp goes over the edge of the end of the chromosome
awk -v OFS="\t" '{ if ($3 == "mRNA") print $1, $4, $5 }' \
../../wbnu.maker.output/wbnu_rnd1.all.maker.noseq.gff | \
awk -v OFS="\t" '{ if ($2 < 1000) print $1, "0", $3+1000; else print $1, $2-1000, $3+1000 }' | \
bedtools getfasta -fi /home/jmanthey/references/wbnu.fasta \
-bed - -fo wbnu_rnd1.all.maker.transcripts1000.fasta


# always activate before running (other conda environments have versions of busco that are broken)
source activate busco

# before running busco edit and export the path to the config file
export BUSCO_CONFIG_FILE="/home/jmanthey/busco/config.ini"

run_busco --in wbnu_rnd1.all.maker.transcripts1000.fasta \
--out busco_output_round1 --lineage_path /home/jmanthey/busco/tetrapoda_odb9/ \
--mode genome -sp human -c 126 -z --long --augustus_parameters='--progress=true'

# job timed out, restarting on xlquanah

#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=augustus
#SBATCH --partition=quanah
#SBATCH --nodes=1
#SBATCH --ntasks=36
#SBATCH --time=120:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --qos=xlquanah

source activate busco

run_busco --in wbnu_rnd1.all.maker.transcripts1000.fasta \
--out busco_output_round1 --lineage_path /home/jmanthey/busco/tetrapoda_odb9/ \
--mode genome -sp human -c 36 -z --long --augustus_parameters='--progress=true' --restart


# when the busco run finishes (this took multiple restarts and ~5 days on 36 processors)
cd /lustre/scratch/jmanthey/00_wbnu_maker/augustus/round1/run_busco_output_round1/augustus_output/retraining_parameters
rename 'BUSCO_busco_output_round1_3344484509' 'Sitta' *
sed -i 's/BUSCO_busco_output_round1_3344484509/Sitta/g' Sitta_parameters.cfg
sed -i 's/BUSCO_busco_output_round1_3344484509/Sitta/g' Sitta_parameters.cfg.orig1 
mkdir $AUGUSTUS_CONFIG_PATH/species/Sitta
cp Sitta_* $AUGUSTUS_CONFIG_PATH/species/Sitta/




