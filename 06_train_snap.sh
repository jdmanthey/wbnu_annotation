# some parts of code from https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2

mkdir -p snap/round1
cd snap/round1

# try various filters and rename and collect stats 
maker2zff -n -d ../../wbnu.maker.output/wbnu_master_datastore_index.log
mv genome.ann wbnu_rnd1.zff.no_filter.ann 
mv genome.dna wbnu_rnd1.zff.no_filter.dna
fathom wbnu_rnd1.zff.no_filter.ann wbnu_rnd1.zff.no_filter.dna -gene-stats > gene-stats_no_filter.log 2>&1

maker2zff -c 0 -e 0 -d ../../wbnu.maker.output/wbnu_master_datastore_index.log
mv genome.ann wbnu_rnd1.zff.aed0.5.ann 
mv genome.dna wbnu_rnd1.zff.aed0.5.dna
fathom wbnu_rnd1.zff.aed0.5.ann wbnu_rnd1.zff.aed0.5.dna -gene-stats > gene-stats_aed0.5.log 2>&1

maker2zff -x 0.25 -l 50 -c 0 -e 0 -d ../../wbnu.maker.output/wbnu_master_datastore_index.log
mv genome.ann wbnu_rnd1.zff.aed0.25_length50.ann 
mv genome.dna wbnu_rnd1.zff.aed0.25_length50.dna
fathom wbnu_rnd1.zff.aed0.25_length50.ann wbnu_rnd1.zff.aed0.25_length50.dna -gene-stats > gene-stats_aed0.25_length50.log 2>&1

maker2zff -x 0.2 -l 75 -c 0 -e 0 -d ../../wbnu.maker.output/wbnu_master_datastore_index.log
mv genome.ann wbnu_rnd1.zff.aed0.2_length75.ann 
mv genome.dna wbnu_rnd1.zff.aed0.2_length75.dna
fathom wbnu_rnd1.zff.aed0.2_length75.ann wbnu_rnd1.zff.aed0.2_length75.dna -gene-stats > gene-stats_aed0.2_length75.log 2>&1

maker2zff -x 0.1 -l 100 -c 0 -e 0 -d ../../wbnu.maker.output/wbnu_master_datastore_index.log
mv genome.ann wbnu_rnd1.zff.aed0.1_length100.ann 
mv genome.dna wbnu_rnd1.zff.aed0.1_length100.dna
fathom wbnu_rnd1.zff.aed0.1_length100.ann wbnu_rnd1.zff.aed0.1_length100.dna -gene-stats > gene-stats_aed0.1_length100.log 2>&1


# check the stats of each and choose the 2nd to last one because it reduced the number of gene models to a reasonable number ~ 9500
fathom wbnu_rnd1.zff.aed0.2_length75.ann wbnu_rnd1.zff.aed0.2_length75.dna -validate > validate.log 2>&1

# categorize with the two files and the surrounding 1000 bp each side
fathom wbnu_rnd1.zff.aed0.2_length75.ann wbnu_rnd1.zff.aed0.2_length75.dna -categorize 1000 > categorize.log 2>&1
fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log 2>&1

# create the training parameters
mkdir params
cd params
forge ../export.ann ../export.dna > ../forge.log 2>&1
cd ..

# make the HMM
hmm-assembler.pl wbnu_rnd1.zff.aed0.2_length75 params > wbnu_rnd1.zff.aed0.2_length75.hmm

