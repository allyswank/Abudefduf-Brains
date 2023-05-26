#!/bin/bash

module load blast+/2.11.0_gcc9

makeblastdb -in uniprot_sprot.fasta -input_type fasta -dbtype prot -out Uniprot.fasta.DB

blastx -db Uniprot.fasta.DB -query Clown_unannotated_contigs.fa -out uniprot_sprot.OUT -evalue 1e-10 -num_threads 16 -outfmt '6 std qcovhsp' -max_target_seqs 5

