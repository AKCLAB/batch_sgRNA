README
=======
WORKFLOW TO IDENTIFY ncRNA

list of subworkflow
# bash input_primer.sh -bedfile ncrna_teste.bed -reffasta LBRAZ_M2903.Dec2022.fasta -size_chr 100 -size_target 50 -outname targetprimer
# bash input_primer.sh -bedfile ncrna_teste.bed -reffasta LBRAZ_M2903.Dec2022.fasta
1. input_primer.sh : 

Pipeline Overview:

Reads the input BED file.
Generates a new BED file with extended coordinates.
Extracts FASTA sequences using a default extension of 150 bp.
Indexes the reference FASTA genome.
Predicts sgRNAs for each target sequence in both the 5′ upstream and 3′ downstream regions of the gene.
Compiles a table with the top 3 sgRNAs predicted by CCTop for each target gene.

Parameters
- size_chr 100
- size_target = 50 (default). This value can be adjusted if necessary: use 50 when the knockout should in some cases occur within the target gene; use 0 if the knockout does not need to fall inside the gene.

Tools: 
cctop 
bowtie2


Inputs: 
- genome reference 
- BED file

Test

bash batch_sgrna.sh -bedfile /cctop_standalone/batch_sgRNA/test/ncrna_teste.bed -reffasta /cctop_standalone/batch_sgRNA/test/LBRAZ_M2903.Dec2022.fasta
genome ref: LBRAZ_M2903.Dec2022.fasta
bed file: ncrna_teste.bed

Outputs:
- outname_initial.bed
- outname_final.bed
- outname_initial.fasta
- outname_final.fasta
- output_cctop/*.xls
- output_cctop/*.bed
- output_cctop/*.fasta
- output_cctop/cctop_listtop3.tsv

python3 posprocessing_outcctop.py

#fiter >900 Candidates are scored from 1000 - suggested best choice to 0 - worst choice. This score takes into account the number of off-targets in the genome, their quality, i.e. number of mismatches and position with respect to the PAM, and the distance to gene exons. 
