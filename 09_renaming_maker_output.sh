# create naming table
maker_map_ids --prefix wbnu --justify 5  wbnu_round2.all.maker.gff > wbnu_round2.all.maker.name.map

# replace names in GFF files
map_gff_ids wbnu_round2.all.maker.name.map wbnu_round2.all.maker.gff
map_gff_ids wbnu_round2.all.maker.name.map wbnu_round2.all.maker.noseqs.gff

# replace names in FASTA headers
map_fasta_ids wbnu_round2.all.maker.name.map wbnu.all.maker.transcripts.fasta
map_fasta_ids wbnu_round2.all.maker.name.map wbnu.all.maker.proteins.fasta

# extract the cds for the entire gff
cat wbnu_round2.all.maker.noseqs.gff | awk '{ if ($3 == "CDS") print $0 }' > wbnu_round2_cds.gff

# extract sequences (make sure to force strandedness)
bedtools getfasta -s -fi /home/jmanthey/references/wbnu.fasta -bed wbnu_round2_cds.gff > wbnu_cds.fasta
