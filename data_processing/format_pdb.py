#write a standard pdb output to be used for phenix refinement
from Bio.PDB import MMCIFParser
from Bio.PDB import PDBParser
from collections import  defaultdict
import os
#follow standard order
_base_atoms_rna = ["P", "C5'","O5'",  "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'"]
_base_atoms_dna = ["P", "C5'","O5'",  "C4'", "O4'", "C3'", "O3'", "C2'", "C1'"]
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
}
_base_atoms_end = ["OP1","OP2","OP3"]
map_dict={"P":"P","C5'":"C","O5'":"O","C4'":"C","O4'":"O","C3'":"C","O3'":"O","C2'":"C","O2'":"O","C1'":"C",
          "N1":"N","C2":"C","O2":"O","N3":"N","C4":"C","N4":"N","C5":"C","C6":"C","O4":"O","N9":"N","C8":"C",
          "N7":"N","N6":"N","O6":"O","OP1":"O","OP2":"O1-","OP3":"O","C7":"C","N2":"N"}
def format_pdb(input_pdb_path,output_pdb_path):
    parser = PDBParser()
    structure = parser.get_structure("input",input_pdb_path)
    with open(output_pdb_path,"w") as wfile:
        atomid=1
        nucid=1
        for model in structure.get_list():
            for chain in model.get_list():
                chain_id = chain.get_id()

                if nucid>9999:
                    nucid=9999
                for residue in chain.get_list():
                    res_name=residue.get_resname().replace(" ","")
                    current_atom_info={} #[atom_name]:[information]
                    for atom in residue.get_list():
                        atom_name = atom.get_fullname().replace(" ","")
                        atom_coord = atom.get_coord()
                        format_coord = []
                        atom_coord = str(atom_coord)
                        atom_coord=atom_coord.replace("[","")
                        atom_coord=atom_coord.replace("]","")
                        atom_coord_split = atom_coord.split()
                        for k in range(3):
                            format_coord.append(float(atom_coord_split[k]))
                        current_atom_info[atom_name]=format_coord
                    if "D" in res_name:
                        _base_atoms_order = _base_atoms_dna
                    else:
                        _base_atoms_order = _base_atoms_rna
                    for atom_name in _base_atoms_order:
                        if atomid>99999:
                            atomid=99999
                        if atom_name in current_atom_info:
                            format_coord=current_atom_info[atom_name]
                            line=""
                            line += "ATOM%7d %-4s %3s%2s%4d    " % (atomid, atom_name,res_name, chain_id,nucid)
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
                            line=""
                            line += "ATOM%7d %-4s %3s%2s%4d    " % (atomid, atom_name,res_name, chain_id,nucid)
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
                            line=""
                            line += "ATOM%7d %-4s %3s%2s%4d    " % (atomid, atom_name,res_name, chain_id,nucid)
                            line = line + "%8.3f%8.3f%8.3f%6.2f%6.2f" % (format_coord[0],format_coord[1],format_coord[2], 1.0, 0)
                            line+=" "*11
                            line+=map_dict[atom_name]+"\n"
                            wfile.write(line)
                            atomid+=1
                    nucid+=1
    return output_pdb_path


if __name__ == "__main__":
    import sys
    format_pdb(sys.argv[1],sys.argv[2])




