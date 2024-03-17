

def Extract_CIF_coord(input_cif_path,filter_list=None):
    # cif_name = os.path.split(input_cif_path)[1][:-4]
    # parser = MMCIFParser()
    # structure = parser.get_structure(cif_name, input_cif_path)
    # for model in structure.get_list():
    #     for chain in model.get_list():
    #         for residue in chain.get_list():
    #read by our own parser


    Coordinate_List=[]
    map_dict={"A":0,"U":1,"T":1,"C":2,"G":3,"DA":0,"DU":1,"DT":1,"DC":2,"DG":3}

    with open(input_cif_path,'r') as rfile:
        for line in rfile:
            if line.startswith("_entry.id"):

                current_id = line.split()[1]


            if line.startswith("ATOM"):
                split_result = line.split()
                coordinates = [float(split_result[10]),
                        float(split_result[11]),
                            float(split_result[12]),]
                atom_name = split_result[3]
                # nuc_type = split_result[5]
                # chain_id = split_result[6]
                # nuc_id = int(split_result[8])
                # pred_label = map_dict[nuc_type]
                # current_score= float(split_result[14])
                if filter_list is not None and atom_name not in filter_list:
                    continue
                Coordinate_List.append(coordinates)

    return Coordinate_List

def cif2pdb(input_cif_path,final_pdb_path):
    begin_check=False
    block_list=[]
    with open(input_cif_path,'r') as rfile:
        for line in rfile:
            if "loop_" in line:
                begin_check=True
                continue

            if begin_check and "_atom_site" in line:
                block_list.append(line.strip("\n").replace(" ",""))
                continue
            if begin_check and "_atom_site" not in line:
                begin_check=False
    atom_ids = block_list.index('_atom_site.id')
    try:
        seq_ids = block_list.index('_atom_site.label_seq_id')
    except:
        seq_ids = block_list.index('_atom_site.auth_seq_id')
    try:
        chain_ids = block_list.index("_atom_site.auth_asym_id")
    except:
        chain_ids = block_list.index("_atom_site.label_asym_id")
    atom_type_ids = block_list.index("_atom_site.label_atom_id")
    res_name_ids = block_list.index("_atom_site.label_comp_id")
    x_ids = block_list.index("_atom_site.Cartn_x")
    y_ids = block_list.index("_atom_site.Cartn_y")
    z_ids = block_list.index("_atom_site.Cartn_z")
    tmp_chain_list=[chr(i) for i in range(ord('A'), ord('Z') + 1)]  # uppercase letters
    tmp_chain_list.extend([chr(i) for i in range(ord('a'), ord('z') + 1)])  # lowercase letters
    #pre filter chain name list
    with open(input_cif_path,'r') as rfile:
        for line in rfile:
            if len(line)>4 and line[:4]=="ATOM":
                split_info=line.strip("\n").split()
                current_chain = split_info[chain_ids]
                if current_chain in tmp_chain_list:
                    tmp_chain_list.remove(current_chain)
    print("remain waiting assign chain number %d"%len(tmp_chain_list))
    chain_map_dict = {}
    with open(input_cif_path,'r') as rfile:
        with open(final_pdb_path,'w') as wfile:

            for line in rfile:
                if len(line)>4 and line[:4]=="ATOM":
                    split_info=line.strip("\n").split()
                    current_chain = split_info[chain_ids]
                    current_atom_index = int(split_info[atom_ids])
                    current_atom_name = split_info[atom_type_ids]
                    current_res_index = int(split_info[seq_ids])
                    current_res_name = split_info[res_name_ids]
                    current_x = float(split_info[x_ids])
                    current_y = float(split_info[y_ids])
                    current_z = float(split_info[z_ids])
                    if len(current_chain)>1: #replace with a temporary id
                        if current_chain in chain_map_dict:
                            current_chain = chain_map_dict[current_chain]
                        else:
                            remain_select_list=[x for x in tmp_chain_list if x not in list(chain_map_dict.values())]
                            chain_map_dict[current_chain]=remain_select_list[0]
                            current_chain = chain_map_dict[current_chain]
                    if current_res_index>9999:
                        current_res_index=9999
                    if current_atom_index>9999999:
                        current_atom_index=9999999
                    wline=""
                    wline += "ATOM%7d %-4s %3s%2s%4d    " % (current_atom_index, current_atom_name,
                                                             current_res_name, current_chain,current_res_index)
                    wline = wline + "%8.3f%8.3f%8.3f%6.2f\n" % (current_x,current_y,current_z, 1.0)
                    wfile.write(wline)