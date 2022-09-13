
from collections import  defaultdict
def read_fasta(input_fasta_path):
    #format should be
    #>chain_id
    #sequence
    dna_rna_set = {"A":0, "U":1, "T":1, "C":2, "G":3}
    chain_dict=defaultdict(list)#key: chain, value: nuc sequence
    with open(input_fasta_path,'r') as file:
        for line in file:
            if line[0]==">":
                current_id = line.strip("\n")
                current_id = current_id.replace(">","")
                if "|" in current_id:
                    current_id = current_id.split("|")[0]
                if "_" in current_id:
                    current_id = current_id.split("_")[1]

            else:
                line=line.strip("\n")
                for item in line:
                    chain_dict[current_id].append(dna_rna_set[item])
    return chain_dict

def read_dna_label(chain_dict):
    dna_label=False
    for key in chain_dict:
        seq = chain_dict[key]
        for item in key:
            if item=="T":
                dna_label=True
    return dna_label


