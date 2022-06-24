

def Show_Graph_Connect(coord_list,edge_list,save_path):
    Natm =1
    with open(save_path,'w') as file:
        file.write('MODEL\n')
        for i in range(len(coord_list)):
            line = ''
            tmp=coord_list[i]
            tmp_chain='A'
            line += "ATOM%7d  %3s %3s%2s%4d    " % (Natm, "CA ", "ALA", " " + tmp_chain, 1)
            line += "%8.3f%8.3f%8.3f%6.2f%6.8f\n" % (
             tmp[0], tmp[1], tmp[2], 1.0, 1.0)
            Natm += 1
            file.write(line)
        for i in range(len(edge_list)):
            nid1,nid2=edge_list[i]
            line = "BOND %d %d\n" % (nid1,nid2)
            file.write(line)

def Show_Bfactor_cif(name,coord_list,save_path,base_prob):
    Natm =1
    chain_dict ={0:"A",1:"B",2:"C",3:"D",4:"E",5:"F",6:"G",7:"H"}
    with open(save_path,'w') as file:
        line = 'data_%s\n'%name
        line += "#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n" \
                   "_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n"\
                    "_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n"\
            "_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n"\
            "_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n_atom_site.auth_seq_id\n_atom_site.auth_asym_id\n"\
            "_atom_site.pdbx_PDB_model_num\n"
        file.write(line)
        for i in range(len(coord_list)):
            tmp=coord_list[i]
            current_prob = base_prob[i]
            tmp_chain="A"
            line =""
            line += "ATOM %-10d C %-3s . %-3s %-2s . %-10d  .  " % (Natm, "CA ", "ALA", " " + tmp_chain, Natm)
            line += "%-8.3f %-8.3f %-8.3f %-6.2f %-6.8f %-10d %-2s 1 \n" % (
             tmp[0], tmp[1], tmp[2], 1.0, current_prob,Natm,tmp_chain)
            Natm += 1
            file.write(line)
