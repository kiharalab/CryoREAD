#write a standard pdb output to be used for phenix refinement
from Bio.PDB import MMCIFParser
from Bio.PDB import PDBParser
from collections import  defaultdict
import os
#follow standard order
_base_atoms_rna = ["P", "C5'","O5'",  "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'"]
_base_atoms_dna = ["P", "C5'","O5'",  "C4'", "O4'", "C3'", "O3'", "C2'", "C1'"]
_base_atoms_protein = ['N',"CA","C","O"]
residue_atoms_order = {
    'A': ["N1", "C2", "N3", "C4","C5","C6", "N6",  "N7", "C8","N9"  ],
    'G':  ["N1", "C2", "N2", "N3", "C4","C5", "C6", "O6","N7", "C8","N9"],
    'C': ["N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"],
    'U': ["N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"],
    'T': ["N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6",'C7'],
    'DA': ["N1", "C2", "N3", "C4","C5","C6", "N6",  "N7", "C8","N9"],
    'DG':  ["N1", "C2", "N2", "N3", "C4","C5", "C6", "O6","N7", "C8","N9"],
    'DC': ["N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"],
    'DT':["N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6",'C7'],
    "ALA": ["CB"],
    "ARG": ["CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"],
    "ASN": ["CB", "CG", "OD1", "ND2"],
    "ASP": ["CB", "CG", "OD1", "OD2"],
    "CYS": ["CB", "SG"],
    "GLU": ["CB", "CG", "CD", "OE1", "OE2"],
    "GLN": ["CB", "CG", "CD", "OE1", "NE2"],
    "GLY": [],
    "HIS": ["CB", "CG", "ND1", "CD2", "CE1", "NE2"],
    "ILE": ["CB", "CG1", "CG2", "CD1"],
    "LEU": ["CB", "CG", "CD1", "CD2"],
    "LYS": ["CB", "CG", "CD", "CE", "NZ"],
    "MET": ["CB", "CG", "SD", "CE"],
    "PHE": ["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "PRO": ["CB", "CG", "CD"],
    "SER": ["CB", "OG"],
    "THR": ["CB", "OG1", "CG2"],
    "TRP": ["CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
    "TYR": ["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"],
    "VAL": ["CB", "CG1", "CG2"]
}
_base_atoms_end = ["OP1","OP2","OP3"]
map_dict={"P":"P","C5'":"C","O5'":"O","C4'":"C","O4'":"O","C3'":"C","O3'":"O","C2'":"C","O2'":"O","C1'":"C",
          "N1":"N","C2":"C","O2":"O","N3":"N","C4":"C","N4":"N","C5":"C","C6":"C","O4":"O","N9":"N","C8":"C",
          "N7":"N","N6":"N","O6":"O","OP1":"O","OP2":"O1-","OP3":"O","C7":"C","N2":"N","N":"N","CA":"C","O":"O","C":"C"}
for key in residue_atoms_order:
    tmp_list=residue_atoms_order[key]
    for atom_name in tmp_list:
        map_dict[atom_name]=atom_name[0]

def write_res_info(wfile,current_atom_info,atomid,nucid):
    if nucid>9999:
        nucid=9999
    if "C4'" in current_atom_info:
        res_name = current_atom_info["C4'"][3]
    elif "CA" in current_atom_info:
        res_name = current_atom_info["CA"][3]
    else:
        return
    if res_name not in residue_atoms_order:
        print("unrecognized info: ",current_atom_info)
        return
    # res_name=residue.get_resname().replace(" ","")
    # current_atom_info={} #[atom_name]:[information]
    # for atom in residue.get_list():
    #     atom_name = atom.get_fullname().replace(" ","")
    #     atom_coord = atom.get_coord()
    #     format_coord = []
    #     atom_coord = str(atom_coord)
    #     atom_coord=atom_coord.replace("[","")
    #     atom_coord=atom_coord.replace("]","")
    #     atom_coord_split = atom_coord.split()
    #     for k in range(3):
    #         format_coord.append(float(atom_coord_split[k]))
    #     current_atom_info[atom_name]=format_coord
    if "CA" in current_atom_info:
        _base_atoms_order = _base_atoms_protein
    elif "D" in res_name:
        _base_atoms_order = _base_atoms_dna
    else:
        _base_atoms_order = _base_atoms_rna
    for atom_name in _base_atoms_order:
        if atomid>99999:
            atomid=99999
        if atom_name in current_atom_info:
            format_coord=current_atom_info[atom_name]
            assert format_coord[3]==res_name
            line=""
            line += "ATOM%7d %-4s %3s%2s%4d    " % (atomid, atom_name,res_name, format_coord[4],nucid)
            line = line + "%8.3f%8.3f%8.3f%6.2f%6.2f" % (format_coord[0],format_coord[1],format_coord[2], 1.0, 0)
            line+=" "*11
            line+=map_dict[atom_name]+"\n"
            wfile.write(line)
            atomid+=1
    current_base_list = residue_atoms_order[res_name]
    for atom_name in current_base_list:
        if atomid>99999:
            atomid=99999
        if atom_name in current_atom_info:
            format_coord=current_atom_info[atom_name]
            assert format_coord[3]==res_name
            line=""
            line += "ATOM%7d %-4s %3s%2s%4d    " % (atomid, atom_name,res_name, format_coord[4],nucid)
            line = line + "%8.3f%8.3f%8.3f%6.2f%6.2f" % (format_coord[0],format_coord[1],format_coord[2], 1.0, 0)
            line+=" "*11
            line+=map_dict[atom_name]+"\n"
            wfile.write(line)
            atomid+=1
    for atom_name in _base_atoms_end:
        if atomid>99999:
            atomid=99999
        if atom_name in current_atom_info:
            format_coord=current_atom_info[atom_name]
            assert format_coord[3]==res_name
            line=""
            line += "ATOM%7d %-4s %3s%2s%4d    " % (atomid, atom_name,res_name, format_coord[4],nucid)
            line = line + "%8.3f%8.3f%8.3f%6.2f%6.2f" % (format_coord[0],format_coord[1],format_coord[2], 1.0, 0)
            line+=" "*11
            line+=map_dict[atom_name]+"\n"
            wfile.write(line)
            atomid+=1
    nucid+=1
    return atomid,nucid
def format_pdb(input_pdb_path,output_pdb_path,DNA_label=False):
    rna2dna_dict={"A":"DA","U":"DT","C":"DC","G":"DG","T":"DT",
                  "DT":"DT","DA":"DA","DC":"DC","DG":"DG"}
    #parser = PDBParser()
    #structure = parser.get_structure("input",input_pdb_path)
    with open(output_pdb_path,"w") as wfile:
        atomid=1
        nucid=1
        prev_resid=0
        keep_info_dict={}
        with open(input_pdb_path,'r') as rfile:
            for read_line in rfile:
                if (read_line.startswith('ATOM')):
                    chain_name = read_line[21]
                    atom_name = read_line[12:16].replace(" ","")
                    x=float(read_line[30:38])
                    y=float(read_line[38:46])
                    z=float(read_line[46:55])
                    resi=int(read_line[22:26])
                    res_name = read_line[17:20].replace(" ","")
                    if DNA_label:
                        res_name = rna2dna_dict[res_name]
                    if resi!=prev_resid:
                        if len(keep_info_dict)>0:
                            atomid,nucid =write_res_info(wfile,keep_info_dict,atomid,nucid)
                        keep_info_dict={}
                        keep_info_dict[atom_name]=[x,y,z,res_name,chain_name]
                        prev_resid=resi
                    else:
                        keep_info_dict[atom_name]=[x,y,z,res_name,chain_name]
        #for model in structure.get_list():
        #    for chain in model.get_list():
                #chain_id = chain.get_id()


    return output_pdb_path

def remove_op3_pdb(input_pdb_path,output_pdb_path):
    with open(input_pdb_path,'r') as rfile:
        with open(output_pdb_path,'w') as wfile:
            for line in rfile:
                if len(line)>4 and line[:4]=="ATOM":
                    atom_name = line[12:16].replace(" ","")
                    if atom_name!="OP3":
                        wfile.write(line)
                else:
                    wfile.write(line)

if __name__ == "__main__":
    import sys
    format_pdb(sys.argv[1],sys.argv[2])




