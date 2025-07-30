import re
import sys
import pandas as pd
import os

#python process_outcctop.py
def process_file(input_xls):
    global blocks, block_atual
    blocks = []
    block_atual = {}

    filename = os.path.basename(input_xls)
    match = re.match(r'^(.*?(?:lncRNA\d+|ncRNA\d+|LBRM\d{5}|LmjF\d{5})).*?(\d{4,}-\d{4,})\.xls$', filename)
    
    # .*? → Captura o início do nome da amostra de forma preguiçosa.
    # (?:lncRNA\d+|ncRNA\d+|LbrM\d{5}|LmjF\d{5}) → Ponto de corte do nome (último identificador útil).
    # .*?(\d{4,6}-\d{4,6}) → Captura a posição final (ex: 104447-104597).
    # \.xls$	Termina com .xls
    if match:
        sample = match.group(1) # Save the name nome, example LbrM2903_01_v2_pilon_lncRNA21, regex (?:lncRNA\d+|ncRNA\d+|LbrM\d{5}|LmjF\d{5}) 
      # Segundo padrão: bctg000000062629-2494bctg000000062579-2729.xls
    else:
        match = re.match(r'^(bctg\d+-\d+)(bctg\d+-\d+)\.xls$', filename)
        if match:
            sample = match.group(1)  # ou group(2), dependendo do que você quiser extrair
        else:
            sample = filename

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
                                            
            elif line.startswith("Chromosome"): # Skip when to identify the line that start with Chromosome
                info_genomic = True
                continue  # skip header with name Chromosome
            elif re.match(r'^(LbrM|chr|contig|bctg|CM)', line):  #If recognize name chromossome with same contain then save 
                if info_genomic:
                    parts = line.split("\t")
                    block_atual["chr"] = parts[0]
                    block_atual["start"] = parts[1]
                    block_atual["end"] = parts[2]
                    block_atual["strand"] = parts[3]
                    block_atual["target_seq"] = parts[5]
                    block_atual["PAM"] = parts[6]
                    info_genomic = False

        # Salva o último bloco
        if block_atual:
            blocks.append(block_atual)

    df = pd.DataFrame(blocks)

    df["efficiency"] = pd.to_numeric(df["efficiency"], errors='coerce').fillna(0).astype(int)
    df["efficiency_CRISPRater"] = pd.to_numeric(df["efficiency_CRISPRater"], errors='coerce')
    df2 = df[df["efficiency"] > 900] # We had consideration this thershold comparing results with high score using Eukaryotic Pathogen CRISPR guide RNA/DNA Design Tool, comparing with CCTOP the same  top sequences present a efficiency with this minimal threshold
    df_sorted = df2.sort_values(by=['efficiency_CRISPRater', 'efficiency'], ascending=False)
    #top3 = df_sorted.head(3)
    return df_sorted.head(3) # Save in the top 3 lines, df_sorted 

if __name__ == "__main__":
    todos_top3 = []

    for file in os.listdir("."):
        if file.endswith(".xls"):
            top3 = process_file(file)
            if not top3.empty:
                todos_top3.append(top3)

    if todos_top3:
        df_final = pd.concat(todos_top3)
        df_final = df_final.sort_values(by=['name',  'efficiency_CRISPRater' ], ascending=[True,False])
        colunas_em_ordem = ['chr','name','target_seq', 'start','end','strand','PAM','id','efficiency','efficiency_CRISPRater']
        df_final.to_csv("cctop_listtop3.tsv", sep="\t", index=False, columns=colunas_em_ordem)
        print("We have the best sgRNAs: cctop_listtop3.tsv")
    else:
        print("No valid top3 found")
