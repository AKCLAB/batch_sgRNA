# Function to display usage information
usage() {
    echo "Usage: $0 -bedfile <file> -reffasta <file> [-dist_chr <value>] [-dist_target <value>] [-outname <file>]"

    echo "  -bedfile: bed file" # model line file: LbrM2903_01_v2_pilon    104547  104783  LbrM2903_01_v2_pilon_lncRNA21   .       +
    echo "  -reffasta: Reference genome FASTA file"
    echo "  -dist_chr: extended size on the chromosome position (optional, default=100)" # To predict sgRNAs, we need to extend the nucleotide information from the start and end position of a target gene, usually it should match at a maximum distance of 100 nucleotides from the start and end position of the gene.
    echo "  -dist_target: extended size on the target position (optional, default=50)" 
    echo "  -outname: name of output bed (optional, default name=targetprimer)"
    exit 1
}
#bash batch_sgrna.sh -bedfile ncrna_teste.bed -reffasta LBRAZ_M2903.Dec2022.fasta -dist_chr 100 -dist_target 50 -outname targetprimer
# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -bedfile) bedfile="$2"; shift; shift ;;
        -reffasta) fasta="$2"; shift; shift ;;
        -dist_chr) dist_chr="$2"; shift; shift ;;
        -dist_target) dist_target="$2"; shift; shift ;;
        -outname) outname="$2"; shift; shift ;;
        *) usage ;;
    esac
done


# Check if all variables are set
if [ -z "$bedfile" ] || [ -z "$fasta" ]; then
    usage
fi
# Set default values
dist_chr="${dist_chr:-100}"
dist_target="${dist_target:-50}"
outname="${outname:-targetprimer}"

# Now use these variables in your script
echo "bed file: $bedfile"
echo "extended size on the chromosome position: $dist_chr"
echo "extended size on the target position: $dist_target"
echo "Reference FASTA: $fasta"
echo "output bed: $outname"

echo "Extract fasta of initial chomossomal position of ncRNA"
awk -v dist_chr="$dist_chr" -v dist_target="$dist_target" 'BEGIN{OFS="\t"} {
    initial_start = $2 - dist_chr;
    initial_end = $2 + dist_target;
    print $1, initial_start, initial_end, $4, $5, $6
}' "$bedfile" > "${outname}_initial.bed"

echo "Extract fasta of end position of ncRNA"
awk -v dist_chr="$dist_chr" -v dist_target="$dist_target" 'BEGIN{OFS="\t"} {
    final_start = $3 - dist_target;
    final_end = $3 + dist_chr;
    print $1, final_start, final_end, $4, $5, $6
}' "$bedfile" > "${outname}_final.bed"

echo "Extraction of fasta from new position according to target"
bedtools getfasta -fi "$fasta" -bed "${outname}_initial.bed" -fo "${outname}_initial.fasta" -name+
bedtools getfasta -fi "$fasta" -bed "${outname}_final.bed" -fo "${outname}_final.fasta" -name+

# verify the index files or run reference index 
ref_name=$(basename "$fasta" .fasta)

# Verify the index file
if ! ls "${ref_name}"*.ebwt 1> /dev/null 2>&1; then
    echo "Running bowtie-build"
    bowtie-build -r -f "$fasta" "$ref_name"
    echo "Index created successfully"
else
    echo "skipping samtools"
fi

path_script="$(dirname "$(realpath "$0")")"
#mkdir output_cctop
# Ensure the output directory exists
if [ ! -d "output_cctop_initial" ]; then
    mkdir "output_cctop_initial"
fi

if [ ! -d "output_cctop_final" ]; then
    mkdir "output_cctop_final"
fi

echo "Running cctop"
cctop --input "${outname}_initial.fasta" --index "$ref_name" --output "${path_script}/output_cctop_initial"
cctop --input "${outname}_final.fasta" --index "$ref_name" --output "${path_script}/output_cctop_final"

echo "pos-processing output cctop for all genes"

python3 "${path_script}/posprocessing_outcctop.py" "${path_script}/output_cctop_initial" "${path_script}/output_cctop_final" "${path_script}/cctop_listtop.tsv"
