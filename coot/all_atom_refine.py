
import os
#work with refine residues


def real_space_atom_refine(input_map,input_pdb,contour):
    output_pdb = input_pdb[:-4]+"_cootrefine.pdb"
    structure = read_pdb(input_pdb)
    map_imol = handle_read_ccp4_map(input_map, 0)
    set_contour_level_absolute(map_imol,contour)
    res_step=1

    residue_list=[]
    for chain_id in chain_ids(structure):
        n_residues = chain_n_residues(chain_id,structure)
        print "There are %(a)i residues in chain %(b)s" % {"a":n_residues,"b":chain_id}

        for serial_number in range(0, n_residues, res_step):
            res_name = resname_from_serial_number(structure, chain_id, serial_number)
            res_no = seqnum_from_serial_number(structure, chain_id, serial_number)
            ins_code = insertion_code_from_serial_number(structure, chain_id, serial_number)
            residue_list.append( [chain_id,res_no,ins_code])
                # active_atom = active_residue()
             # centred_residue = active_atom[1:4]
    with AutoAccept():
        refine_residues(structure,residue_list)
    write_pdb_file(structure, output_pdb)
    close_molecule(structure)
    close_molecule(map_imol)
    coot_real_exit(0)

def real_space_atom_refine_output(input_map,input_pdb,output_pdb,contour):
    #output_pdb = input_pdb[:-4]+"_cootrefine.pdb"
    structure = read_pdb(input_pdb)
    map_imol = handle_read_ccp4_map(input_map, 0)
    set_contour_level_absolute(map_imol,contour)
    res_step=1

    residue_list=[]
    for chain_id in chain_ids(structure):
        n_residues = chain_n_residues(chain_id,structure)
        print "There are %(a)i residues in chain %(b)s" % {"a":n_residues,"b":chain_id}

        for serial_number in range(0, n_residues, res_step):
            res_name = resname_from_serial_number(structure, chain_id, serial_number)
            res_no = seqnum_from_serial_number(structure, chain_id, serial_number)
            ins_code = insertion_code_from_serial_number(structure, chain_id, serial_number)
            residue_list.append( [chain_id,res_no,ins_code])
                # active_atom = active_residue()
             # centred_residue = active_atom[1:4]
    with AutoAccept():
        refine_residues(structure,residue_list)
    write_pdb_file(structure, output_pdb)
    close_molecule(structure)
    close_molecule(map_imol)
    coot_real_exit(0)

def real_space_atom_fragmentrefine_output(input_map,input_pdb,output_pdb,contour):
    #output_pdb = input_pdb[:-4]+"_cootrefine.pdb"
    structure = read_pdb(input_pdb)
    map_imol = handle_read_ccp4_map(input_map, 0)
    set_contour_level_absolute(map_imol,contour)
    res_step=1

    residue_list=[]
    #begin_flag=True
    prev_residue =-1
    prev_chain="GG"
    for chain_id in chain_ids(structure):
        n_residues = chain_n_residues(chain_id,structure)
        print "There are %(a)i residues in chain %(b)s" % {"a":n_residues,"b":chain_id}

        for serial_number in range(0, n_residues, res_step):
            res_name = resname_from_serial_number(structure, chain_id, serial_number)
            res_no = seqnum_from_serial_number(structure, chain_id, serial_number)
            ins_code = insertion_code_from_serial_number(structure, chain_id, serial_number)




            if chain_id==prev_chain and abs(res_no-prev_residue)==1:
                residue_list.append([chain_id,res_no,ins_code])
                prev_residue=res_no
                prev_chain=chain_id
            else:
                if len(residue_list)!=0:
                    print("refine structures: ",residue_list)
                    with AutoAccept():
                        refine_residues(structure,residue_list)
                residue_list=[]
                residue_list.append([chain_id,res_no,ins_code])
                prev_residue=res_no
                prev_chain=chain_id
                # active_atom = active_residue()
             # centred_residue = active_atom[1:4]
        if len(residue_list)!=0:
            print("refine structures: ",residue_list)
            with AutoAccept():
                refine_residues(structure,residue_list)
    write_pdb_file(structure, output_pdb)
    close_molecule(structure)
    close_molecule(map_imol)
    coot_real_exit(0)

def real_space_atom_residuerefine_output(input_map,input_pdb,output_pdb,contour):
    #output_pdb = input_pdb[:-4]+"_cootrefine.pdb"
    structure = read_pdb(input_pdb)
    map_imol = handle_read_ccp4_map(input_map, 0)
    set_contour_level_absolute(map_imol,contour)
    res_step=1

    residue_list=[]
    #begin_flag=True
    prev_residue =-1
    prev_chain="GG"
    for chain_id in chain_ids(structure):
        n_residues = chain_n_residues(chain_id,structure)
        print "There are %(a)i residues in chain %(b)s" % {"a":n_residues,"b":chain_id}

        for serial_number in range(0, n_residues, res_step):
            res_name = resname_from_serial_number(structure, chain_id, serial_number)
            res_no = seqnum_from_serial_number(structure, chain_id, serial_number)
            ins_code = insertion_code_from_serial_number(structure, chain_id, serial_number)


            with AutoAccept():
                refine_residues(structure,[[chain_id,res_no,ins_code]])

    write_pdb_file(structure, output_pdb)
    close_molecule(structure)
    close_molecule(map_imol)
    coot_real_exit(0)

def real_space_atom_3residuerefine_output(input_map,input_pdb,output_pdb,contour):
    #output_pdb = input_pdb[:-4]+"_cootrefine.pdb"
    structure = read_pdb(input_pdb)
    map_imol = handle_read_ccp4_map(input_map, 0)
    set_contour_level_absolute(map_imol,contour)
    res_step=1
    all_list=[]
    residue_list=[]
    #begin_flag=True
    prev_residue =-1
    prev_chain="GG"
    for chain_id in chain_ids(structure):
        n_residues = chain_n_residues(chain_id,structure)
        print "There are %(a)i residues in chain %(b)s" % {"a":n_residues,"b":chain_id}

        for serial_number in range(0, n_residues, res_step):
            res_name = resname_from_serial_number(structure, chain_id, serial_number)
            res_no = seqnum_from_serial_number(structure, chain_id, serial_number)
            ins_code = insertion_code_from_serial_number(structure, chain_id, serial_number)




            if chain_id==prev_chain and abs(res_no-prev_residue)==1:
                residue_list.append([chain_id,res_no,ins_code])
                prev_residue=res_no
                prev_chain=chain_id
            else:
                if len(residue_list)!=0:
                    print("collecting structures: ",residue_list)
                    all_list.append(residue_list)
                residue_list=[]
                residue_list.append([chain_id,res_no,ins_code])
                prev_residue=res_no
                prev_chain=chain_id
                # active_atom = active_residue()
             # centred_residue = active_atom[1:4]
        if len(residue_list)!=0:
            print("collecting structures: ",residue_list)
            all_list.append(residue_list)

    for residue_list in all_list:
        if len(residue_list)<=3:
            print("refine structures: ",residue_list)
            with AutoAccept():
                refine_residues(structure,residue_list)
        else:

            for k in range(0,len(residue_list)-3):
                print("refine structures: ",residue_list[k:k+3])
                with AutoAccept():
                    refine_residues(structure,residue_list[k:k+3])
    write_pdb_file(structure, output_pdb)
    close_molecule(structure)
    close_molecule(map_imol)
    coot_real_exit(0)

def real_space_atom_Bchainrefine_output(input_map,input_pdb,output_pdb,contour):
    #output_pdb = input_pdb[:-4]+"_cootrefine.pdb"
    structure = read_pdb(input_pdb)
    map_imol = handle_read_ccp4_map(input_map, 0)
    set_contour_level_absolute(map_imol,contour)
    res_step=1

    residue_list=[]
    for chain_id in chain_ids(structure):
        if chain_id!='B':
            continue
        n_residues = chain_n_residues(chain_id,structure)
        print "There are %(a)i residues in chain %(b)s" % {"a":n_residues,"b":chain_id}

        for serial_number in range(0, n_residues, res_step):
            res_name = resname_from_serial_number(structure, chain_id, serial_number)
            res_no = seqnum_from_serial_number(structure, chain_id, serial_number)
            ins_code = insertion_code_from_serial_number(structure, chain_id, serial_number)
            residue_list.append( [chain_id,res_no,ins_code])
                # active_atom = active_residue()
             # centred_residue = active_atom[1:4]
    with AutoAccept():
        refine_residues(structure,residue_list)
    write_pdb_file(structure, output_pdb)
    close_molecule(structure)
    close_molecule(map_imol)
    coot_real_exit(0)

def real_space_atom_range_refine_output(input_map,input_pdb,output_pdb,contour,refine_index):
    #output_pdb = input_pdb[:-4]+"_cootrefine.pdb"
    refine_start, refine_end = refine_index
    structure = read_pdb(input_pdb)
    map_imol = handle_read_ccp4_map(input_map, 0)
    set_contour_level_absolute(map_imol,contour)
    res_step=1

    residue_list=[]
    for chain_id in chain_ids(structure):

        n_residues = chain_n_residues(chain_id,structure)
        print "There are %(a)i residues in chain %(b)s" % {"a":n_residues,"b":chain_id}

        for serial_number in range(0, n_residues, res_step):
            res_name = resname_from_serial_number(structure, chain_id, serial_number)
            res_no = seqnum_from_serial_number(structure, chain_id, serial_number)
            ins_code = insertion_code_from_serial_number(structure, chain_id, serial_number)
            if int(res_no)>=refine_start and int(res_no)<=refine_end:
                residue_list.append( [chain_id,res_no,ins_code])
                # active_atom = active_residue()
             # centred_residue = active_atom[1:4]
    print("refine list",residue_list)
    with AutoAccept():
        refine_residues(structure,residue_list)
    write_pdb_file(structure, output_pdb)
    close_molecule(structure)
    close_molecule(map_imol)
    coot_real_exit(0)
