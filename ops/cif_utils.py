

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
