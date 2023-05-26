#!/bin/bash

module load blast+/2.11.0_gcc9

makeblastdb -in Amphiprion_coding_seq.fa -input_type fasta -dbtype nucl -out Amphiprion_percula.fa.DB

blastn -db Amphiprion_percula.fa.DB -query Asax_truncnames.fa -out Nemo_full_gene.OUT -evalue 1e-10 -num_threads 16 -outfmt '6 std qcovhsp' -max_target_seqs 5

