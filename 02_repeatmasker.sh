# use a custom repeatmasker database to annotate the genome for TEs 
# includes the:
# RepBase vertebrate database v24.03 sequences
# certhia americana custom repeatmodeler sequences (doi: 10.1093/gbe/evab120)
# colaptes auratus custom repeatmodeler sequences (doi:10.1093/g3journal/jkaa026)

# run repeat masker v1.332
cd /home/jmanthey/references
RepeatMasker -pa 24 -s -lib ~/RepeatMasker/Libraries/custom_library_certhia_colaptes.fa wbnu_genome_NewNames__final_assembly.fasta
