
from collections import  defaultdict
def read_fasta(input_fasta_path,dna_check=False):
    #format should be
    #>chain_id
    #sequence
    dna_rna_set = {"A":0, "U":1, "T":1, "C":2, "G":3}
    chain_dict=defaultdict(list)#key: chain, value: nuc sequence
    current_id=None

    tmp_chain_list=[chr(i) for i in range(ord('A'), ord('Z') + 1)]  # uppercase letters
    tmp_chain_list.extend([chr(i) for i in range(ord('a'), ord('z') + 1)])  # lowercase letters
    read_chain=False
    dna_label=True
    with open(input_fasta_path,'r') as file:
        for line in file:
            if line[0]==">":
                current_id = line.strip("\n")
                current_id = current_id.replace(">","")
                if "|" in current_id:
                    current_id = current_id.split("|")[0]
                if "_" in current_id:
                    current_id = current_id.split("_")[1]
                read_chain=True
            elif len(line.strip("\n").replace(" ",""))>0:
                if current_id is None or read_chain==False or len(current_id)!=1:
                    visit_set=list(chain_dict.keys())
                    for tmp_chain in tmp_chain_list:
                        if tmp_chain not in visit_set:
                            current_id=tmp_chain
                            break

                line=line.strip("\n").replace(" ","").replace("N","").replace("X","")
                #quick check to see if it includes protein fasta
                dna_seq_flag=True
                count_porition=0
                for item in line:
                    if item in dna_rna_set:
                        count_porition+=1
                count_porition = count_porition/len(line)
                if count_porition<=0.9:
                    dna_seq_flag=False
                if dna_seq_flag:
                    for item in line:
                        if item in dna_rna_set:
                            if item=="U":
                                dna_label=False
                            chain_dict[current_id].append(dna_rna_set[item])
    print("read chain info from fasta:",chain_dict)
    if dna_check:
        return chain_dict,dna_label
    else:
        return chain_dict

# def read_dna_label(chain_dict):
#     dna_label=False
#     for key in chain_dict:
#         seq = chain_dict[key]
#         for item in key:
#             if item=="T":
#                 dna_label=True
#     return dna_label


