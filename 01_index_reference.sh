# move to location of reference
cd references/

#index with samtools and bwa
samtools faidx wbnu.fasta
bwa index wbnu.fasta

# move up a directory to location of picard jar file
cd ../
# make the dict file
java -jar picard.jar CreateSequenceDictionary R=/home/jmanthey/references/wbnu.fasta O=/home/jmanthey/references/wbnu.dict
