

def read_fasta(input_fasta_path):
    #format should be
    #>chain_id
    #sequence
    chain_dict={}#key: chain, value: nuc sequence
    with open(input_fasta_path,'r') as file:
        for line in file:
            if line[0]==">":
                current_id = line.strip("\n")
                current_id = current_id.replace(">","")
            else:
                line=line.strip("\n")
                chain_dict[current_id]=line
    return chain_dict

def read_dna_label(chain_dict):
    dna_label=False
    for key in chain_dict:
        seq = chain_dict[key]
        for item in key:
            if item=="T":
                dna_label=True
    return dna_label


