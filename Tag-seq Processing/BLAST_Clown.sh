#!/bin/bash

module load blast+/2.11.0_gcc9

makeblastdb -in Clown_coding_seq.fa -input_type fasta -dbtype nucl -out Clownfish.fa.DB

blastn -db Clownfish.fa.DB -query Steg_unannotated_contigs.fa -out Clown_full_gene.OUT -evalue 1e-10 -num_threads 16 -outfmt '6 std qcovhsp' -max_target_seqs 5

