#!/bin/sh

indir="/home/aubars001/a-sax-rnaseq/annotations" # directory containing your BLAST output table, just get rid of the text within the quotes and fill in your appropriate directory
blastfile="./uniprot_sprot.OUT" # the name of your output table from BLAST
goinfo="uniprot_GO_terms.csv" # the name of a file containing genes and GO term info from your organism that you BLASTed against.
# The GO info file for this was downloaded by taking a list of uniprot IDs and searching for them on uniprot's retrieve/ID mapping page (uniprot.org/uploadlists). I grabbed the uniprot IDs from the uniprot_sprot_ids.txt table.
# I used awk to parse out the six-character IDs, which I then input into a file to be searched on uniprot. I searched uniprot IDs against uniprot KB, and then on the output page, I went to the "columns" tab to select/deselect what I wanted to include. I added columns for GO term numbers and phylum (to remove any non-vertebrate matches as those seem like they would probably be incorrect). 
# I downloaded this to my local machine as an excel file, did find/replace to remove any commas, and then saved as a CSV to be uploaded back to the ASC.
mincoverage="16" # how much of the "query" sequence you want to be covered by the "match" sequence. This is the value you have to multiply by 3 to get the actual percentage.
minpercentid="60" # how high the percent identity needs to be before you consider something a match
evalue="1e-7" # how low the evalue needs to be before you consider something a match 

# Navigate to the input directory
cd $indir

# This filters your table based on the parameters given above:
# First, it filters by the query coverage, then it takes the output of that filtering and filters for percent identity, finally, it takes the output of that filtering and filters for e-value. 
# The resulting table is written to a file called multiple_matches.txt
cat $blastfile | awk -v mincov="$mincoverage" '$13 > mincov {print $0}' | awk -v minid="$minpercentid" '$3 >= minid {print $0}' | awk -v ev="$evalue" '$11 < ev {print $0}' > uniprot_sprot_multiple_matches.txt

# Next, we will get a list of all the transcripts that have annotation information
# This file will be called annotated_transcripts.txt
cat uniprot_sprot_multiple_matches.txt | awk '{print $1}' | sort | uniq > uniprot_sprot_annotated_transcripts.txt

# Now, for every transcript listed in that file, search for the name of that transcript in your multiple_matches file, and output that into a file called temp.txt
# Sort the matches based on column 12, which is bit-score
# I would rather that this sort by e-value, but sort doesn't work properly with numbers in exponential notation and I don't know how to get around it. Bit score is also a fine metric to sort by though.
# Then, this program takes the top match (based on bit score) and outputs it into a final file called one_annotation_per_transcript
# Remove the temporary files and then move to next iteration of the loop
for transcript in `cat uniprot_sprot_annotated_transcripts.txt`
do
    egrep "$transcript\s" uniprot_sprot_multiple_matches.txt > temp.txt
    sort -nr -k 12 temp.txt > sorted.txt
    head -1 sorted.txt >> uniprot_sprot_one_annotation_per_transcript.txt
    rm temp.txt
    rm sorted.txt
done

# From the one_annotation_per_transcript file, print ONLY the TRINITY_ID, matching SPROT_ID, percent identity, e-value, bit-score, and query coverage
# Replace any spaces in this file with commas, because spaces are stupid and commas work wonders, output this into a file called ids.txt
cat uniprot_sprot_one_annotation_per_transcript.txt | awk '{print $1,$2,$3,$11,$12,$13}' | sed -r 's/\s+/,/g' > uniprot_sprot_ids.txt

# Finally, for each gene in the ids.txt file (contains ensembl and trinity IDs), do the following:
# Search for that gene's SPROT ID in another file, which contains SPROT IDs, gene names, GO terms, and GO accessions
# Output the contents of that search into a temporary file
# Then, for every line in that temporary file, add the associated trinity ID to the beginning of the line, and append this output to annotation_table.csv.
for gene in `cat uniprot_sprot_ids.txt`
do
	TRINITY=`echo "$gene" | awk -F, '{print $1}'`
	PID=`echo "$gene" | awk -F, '{print $3}'`
	EVAL=`echo "$gene" | awk -F, '{print $4}'`
	BIT=`echo "$gene" | awk -F, '{print $5}'`
	COV=`echo "$gene" | awk -F, '{print $6}'`	
	SPROTI=`echo "$gene" | awk -F, '{print $2}' | awk -F'|' '{print $2}'`
	grep "$SPROTI" $goinfo | sed "s/^$SPROTI/$TRINITY,$PID,$EVAL,$BIT,$COV,$SPROTI/g" >> swissprot_annotation_table.csv
done

# The final output will be a file called annotation_table.csv. The columns in this table will be, in the following order:
# Trinity ID, Percent Identity of blast match, e-value of blast match, bit-score of blast match, query coverage of blast match, uniprot id, followed by the matching gene and go term info from your GO term file.

# You can use the following code to find UN-annotated contigs that are not in your final annotation file, so they can then be BLASTed or searched for again with a different gene set (like from another species)
# Fill in "name-of-your-transcriptome" with the name of your fasta file containing your transcriptome assembly
# The annotated transcripts.txt file has been generated in the previous steps

source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load seqkit/0.10.1

seqkit grep -v -f uniprot_sprot_annotated_transcripts.txt Clown_unannotated_contigs.fa -o uniprot_sprot_unannotated_contigs.fa

# The results of this are a file which contains ALL of the matches from SwissProt (NOT JUST THE CHORDATES). 
# This is a problem because BLAST matches to proteins outside of chordata are probably inaccurate. 
# We can filter our results to only include chordates using the following code:

awk -F, '$13 ~ "Chordata" {print}' swissprot_annotation_table.csv > swissprot_annotation_table_chordates.csv

awk -F, '{print $1}' swissprot_annotation_table_chordates.csv > uniprot_sprot_annotated_transcripts_chordates.txt

seqkit grep -v -f uniprot_sprot_annotated_transcripts_chordates.txt Clown_unannotated_contigs.fa -o uniprot_sprot_chordates_unannotated_contigs.fa
