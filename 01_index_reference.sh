# move to location of reference
cd references/

#index with samtools and bwa
samtools faidx wbnu_genome_NewNames__final_assembly.fasta
bwa index wbnu_genome_NewNames__final_assembly.fasta

# move up a directory to location of picard jar file
cd ../
# make the dict file
java -jar picard.jar CreateSequenceDictionary R=/home/jmanthey/references/wbnu_genome_NewNames__final_assembly.fasta O=/home/jmanthey/references/wbnu_genome_NewNames__final_assembly.dict
