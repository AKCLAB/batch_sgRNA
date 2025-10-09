import re
import sys
import pandas as pd
import os
from Bio import SeqIO
#python3 posprocessing_outcctop.py file.fasta output_cctop_initial output_cctop_final cctop_listtop.tsv file.gff
top = 3
def process_file(input_xls, position):
    global blocks, block_atual
    blocks = []
    block_atual = {}

    filename = os.path.basename(input_xls)
    #Extract name from sample file
    #match = re.match(r'^(.*?(?:lncRNA\d+|ncRNA\d+|LBRM\d{5}|LmjF\d{5})).*?(\d{4,}-\d{4,})\.xls$', filename)
    match = re.match(r'^(.*?)(?:_)?((?:lncRNA\d+|ncRNA\d+|LBRM\d{5}|LmjF\d{5})).*?(\d{4,}-\d{4,})\.xls$', filename)
    # .*? → Captura o início do nome da amostra de forma preguiçosa.
    # (?:lncRNA\d+|ncRNA\d+|LbrM\d{5}|LmjF\d{5}) → Ponto de corte do nome (último identificador útil).
    # .*?(\d{4,6}-\d{4,6}) → Captura a posição final (ex: 104447-104597).
    # \.xls$	Termina com .xls 
    if match:
        chromosome = match.group(1)  # ex: LbrM2903_01_v2_pilon ou bctg00000022
        sample = f"{match.group(1)}_{match.group(2)}"  # ex: LbrM2903_01_v2_pilon_lncRNA21
    # Segundo padrão: bctg000000062629-2494bctg000000062579-2729.xls
    else:
        match = re.match(r'^(bctg\d+-\d+)(bctg\d+-\d+)\.xls$', filename)
        if match:
            sample = match.group(1)
        else:
            sample = filename
    info_genomic = False
    with open(input_xls, "r") as f:
        for i, line in enumerate(f):
            if i < 8:
                continue  # skip the first 8 lines
            line = line.strip()
            if line.startswith("T") and "\t" in line: # identify the prediction with start  (ex: T1, T2...)
                if block_atual:  # Save the previous block
                    blocks.append(block_atual)
                    block_atual = {} #save information in library
                parts = line.split("\t")
                block_atual["name"] = sample
                block_atual["id"] = parts[0] # Save the first column of line that start with T1
                #block_atual["sequence"] = parts[1] # Save the second column
                block_atual["efficiency"] = parts[2]
                block_atual["efficiency_CRISPRater"] = parts[-1]
                block_atual["position"] = position
                                            
            elif line.startswith("Chromosome"): # Skip when to identify the line that start with Chromosome
                info_genomic = True
                continue  # skip header with name Chromosome
            if info_genomic and re.match(r'^(?:LbrM|chr|contig|bctg|CM|LmjF|LinJ|Tb)', line):
                parts = line.split("\t")
                # Cheack at least 6 columns
                if len(parts) < 6:
                    #  Ignore
                    continue
                chrom_in_line = parts[0].strip()
                # Continue saving if the chromossome is the same the chromossome from name  file
                if chromosome and chrom_in_line == chromosome:
                    parts = line.split("\t")
                    block_atual["chr"] = parts[0]
                    block_atual["start"] = parts[1]
                    block_atual["end"] = parts[2]
                    block_atual["strand"] = parts[3]
                    block_atual["target_seq"] = parts[5] # sgRNA sequence
                    block_atual["PAM"] = parts[6]
                    info_genomic = False
        # Save the last block
        if block_atual:
            blocks.append(block_atual)
    df = pd.DataFrame(blocks)

    # Load genome into a dictionary to extract 30 nt for homologous recombination following the cut site of the PAM sequence at the 3′ end of the guide sequence
    genome = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    # List to save sequence
    seq1_list = []
    for idx, row in df.iterrows():  
        #Cut site at the PAM sequence on the 3′ end of the guide of downstream target region.
        if row["position"] == "downstream" :
            if row["strand"] == "+" :
                chrom = row["chr"]
                hr1 = int(row["end"]) - 6 
                # extrair fragmentos de 30 nt (ajustando para 0-based do python)
                seq1 = genome[chrom].seq[hr1:hr1+30]
                seq1_list.append(str(seq1))
            elif row["strand"] =="-" :
                chrom = row["chr"]
                hr1 = int(row["start"]) + 5 
                seq1 = genome[chrom].seq[hr1:hr1+30].reverse_complement()
                seq1_list.append(str(seq1))
        #Cut site at the PAM sequence on the 3′ end of the guide of upstream target region.        
        if row["position"] == "upstream" :
            if row["strand"] == "+" :
                chrom = row["chr"]
                hr1 = int(row["end"]) - 6 
                # extrair fragmentos de 30 nt (ajustando para 0-based do python)
                seq1 = genome[chrom].seq[hr1-30:hr1]                
                seq1_list.append(str(seq1))
            elif row["strand"] =="-" :
                chrom = row["chr"]
                hr1 = int(row["start"]) + 5 
                seq1 = genome[chrom].seq[hr1-30:hr1].reverse_complement()
                seq1_list.append(str(seq1))
            # adicionar colunas na tabela
    df["HR"] = seq1_list
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)

    # If GFF is provided, it checks to remove sgRNA overlap CDs or ncRNA.
    if gff and gff.strip():
        print(f"Verificando sobreposição usando {gff} ...")
        # Function to extract the ID from the 9th column in a GFF file
        def extract_info(text):
            id_match = re.search(r'ID=([^;]+)', text)  # extract ID after "ID=" and before ";"
            id_value = id_match.group(1) if id_match else None
            return id_value

        # Creates a Boolean Series of False to indicate global overlap
        df["overlap"] = False
        with open(gff, 'r') as fh_cds:
            # Verify sgRNA overlapping on CDsA
            for cds_row in fh_cds:
                if cds_row.startswith("#"):
                    continue
                # Split the line into columns
                cds_fields = cds_row.strip().split("\t")
                # Check if the line is a CDS or ncRNA
                if cds_fields[2] in ["CDS", "ncRNA"]:
                    cds_chr, cds_coordi, cds_coordf = cds_fields[0],  int(cds_fields[3]), int(cds_fields[4])
                    # if the sgRNA is in the same CDS or ncRNA
                    same_chr = df["chr"] == cds_chr
                    # Verify overlapping
                    mask_overlap = (
                        ((df["start"] >= cds_coordi) & (df["start"] <= cds_coordf)) |
                        ((df["end"] >= cds_coordi) & (df["end"] <= cds_coordf)) |
                        ((cds_coordi >= df["start"]) & (cds_coordi <= df["end"])) |
                        ((cds_coordf >= df["start"]) & (cds_coordf <= df["end"]))
                    )
                    #((df["start"] >= cds_coordi) & (df["start"] <= cds_coordf))→ o start da sgRNA está dentro do CDS (início do sgRNA entre o início e o fim do CDS).
                    #((df["end"] >= cds_coordi) & (df["end"] <= cds_coordf)) → o end do sgRNA está dentro do CDS.
                    #((cds_coordi >= df["start"]) & (cds_coordi <= df["end"]))→ o início do CDS está dentro do intervalo do sgRNA (CDS começa dentro do sgRNA).
                    #((cds_coordf >= df["start"]) & (cds_coordf <= df["end"]))→ o fim do CDS está dentro do intervalo do sgRNA
                    df.loc[same_chr & mask_overlap, "overlap"] = True
        # Keep only those that do NOT overlap with any gene
        df = df[~df["overlap"]].drop(columns="overlap")

    df["efficiency"] = pd.to_numeric(df["efficiency"], errors='coerce').fillna(0).astype(int)
    df["efficiency_CRISPRater"] = pd.to_numeric(df["efficiency_CRISPRater"], errors='coerce')
    df2 = df[df["efficiency"] > 900] # We had consideration this thershold comparing results with high score using Eukaryotic Pathogen CRISPR guide RNA/DNA Design Tool, comparing with CCTOP the same  top sequences present a efficiency with this minimal threshold
    df_sorted = df2.sort_values(by=['efficiency_CRISPRater', 'efficiency'], ascending=False)
    return df_sorted.head(top) # Save in the top 3 lines, df_sorted 

if __name__ == "__main__":
    # Espera 4 ou 5 argumentos
    if len(sys.argv) < 5 or len(sys.argv) > 6:
        print("Uso: python posprocessing_outcctop.py <fasta> <dir_initial> <dir_final> <saida.tsv> <gff>")
        sys.exit(1)
    fasta = sys.argv[1]
    dir_initial = sys.argv[2]
    dir_final = sys.argv[3]
    output_path = sys.argv[4]
    gff = sys.argv[5] if len(sys.argv) == 6 else None  # opcional
    todos_top = []

    # Import initial file
    for file in os.listdir(dir_initial):
        if file.endswith(".xls"):
            top3 = process_file(os.path.join(dir_initial, file), "upstream")
            if not top3.empty:
                todos_top.append(top3)

    # Import final file
    for file in os.listdir(dir_final):
        if file.endswith(".xls"):
            top3 = process_file(os.path.join(dir_final, file), "downstream")
            if not top3.empty:
                todos_top.append(top3)

    if todos_top:
        df_final = pd.concat(todos_top)
        df_final = df_final.sort_values(by=['name', 'efficiency_CRISPRater'], ascending=[True, False])
        colunas_em_ordem = ['chr', 'name', 'target_seq', 'start', 'end', 'strand', 'PAM', 'id', 'efficiency', 'efficiency_CRISPRater', 'position', 'HR']
        df_final.to_csv(output_path, sep="\t", index=False, columns=colunas_em_ordem)
        print("We have the best sgRNAs: {}".format(output_path))

    else:
        print("No valid top found")

