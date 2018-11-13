#! /bin/sh

# Make the blast indexes for the assignment sequences 
makeblastdb -in drosophila.fa -parse_seqids -dbtype prot
makeblastdb -in celegans.fa -parse_seqids -dbtype prot

# Do the reciprocal search of the protein sequences
blastp -query celegans.fa -db drosophila.fa -outfmt 5 -out celegans_drosophila.xml
blastp -query drosophila.fa -db celegans.fa -outfmt 5 -out drosophila_celegans.xml

