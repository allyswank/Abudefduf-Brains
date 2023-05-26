#!/bin/bash

module load blast+/2.11.0_gcc9

makeblastdb -in Chromis_coding_seq.fa -input_type fasta -dbtype nucl -out Spiny_Chromis.fa.DB

blastn -db Spiny_Chromis.fa.DB -query Nemo_unannotated_contigs.fa -out Chromis_full_gene.OUT -evalue 1e-10 -num_threads 16 -outfmt '6 std qcovhsp' -max_target_seqs 5

