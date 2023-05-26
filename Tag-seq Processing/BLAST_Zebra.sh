#!/bin/bash

module load blast+/2.11.0_gcc9

makeblastdb -in Zebra_coding_seq.fa -input_type fasta -dbtype nucl -out Zebrafish.fa.DB

blastn -db Zebrafish.fa.DB -query Chromis_unannotated_contigs.fa -out Zebra_full_gene.OUT -evalue 1e-10 -num_threads 16 -outfmt '6 std qcovhsp' -max_target_seqs 5

