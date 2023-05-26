#!/bin/bash

module load blast+/2.11.0_gcc9

makeblastdb -in Steg_coding_seq.fa -input_type fasta -dbtype nucl -out Stegastes.fa.DB

blastn -db Stegastes.fa.DB -query Zebra_unannotated_contigs.fa -out Steg_full_gene.OUT -evalue 1e-10 -num_threads 16 -outfmt '6 std qcovhsp' -max_target_seqs 5

