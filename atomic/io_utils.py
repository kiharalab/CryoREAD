from collections import defaultdict
import numpy as np
import os
from scipy.spatial.distance import cdist
from atomic.geo_atomic_modeling import build_base_atomic_nuc_rigid,build_base_atomic_nuc_rigid_sugar

def define_base_location_list(frag_info_list,refer_base_location):
    all_frag_coordinate= []
    for i in range(len(frag_info_list)):
        chain_id,current_seq_index,atom_name, cur_atom_position,nuc_type,avg_score =frag_info_list[i]
        if atom_name=="C4'":
            all_frag_coordinate.append(cur_atom_position)
    all_refer_coordinate = []
    all_refer_key = list(refer_base_location.keys())
    for key in all_refer_key:
        split_key = key.split(",")
        cur_coord=[]
        for k in range(3):
            cur_coord.append(float(split_key[k]))
        all_refer_coordinate.append(cur_coord)
    all_refer_coordinate = np.array(all_refer_coordinate)
    refer_distance = cdist(all_frag_coordinate,all_refer_coordinate)
    refer_close = np.argmin(refer_distance,axis=1)
    refer_base_location_list = []
    for k in range(len(all_frag_coordinate)):
        close_key = int(refer_close[k])
        cur_distance = refer_distance[k,close_key]
        assert cur_distance<1#should be the same coord used for
        search_key = all_refer_key[close_key]
        refer_base_location_list.append(refer_base_location[search_key])
    return refer_base_location_list
def Write_Atomic_Fraginfo_cif(name,frag_info_list,refer_base_location,save_path,DNA_Label,map_info_list):
    Natm =1
    pdb_save_path = save_path[:-4]+".pdb"
    with open(pdb_save_path,"w") as file:
        file.write("")
    refer_sugar_base_location = define_base_location_list(frag_info_list,refer_base_location)
    with open(save_path,'w') as file:
        line = 'data_%s\n'%name
        line += "#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n" \
                   "_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n"\
                    "_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n"\
            "_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n"\
            "_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n_atom_site.auth_seq_id\n_atom_site.auth_asym_id\n"\
            "_atom_site.pdbx_PDB_model_num\n"
        file.write(line)
        count_track = 0
        #extract sugar and phosphate_pairs
        Sugar_Refer_Dict={}#use nuc_id to identify Phosphate
        Sugar_Only_Frag_List=[]
        for i in range(len(frag_info_list)):
            chain_id,current_seq_index,atom_name, cur_atom_position,nuc_type,avg_score =frag_info_list[i]
            if atom_name=="C4'":
                Sugar_Only_Frag_List.append(frag_info_list[i])
            else:
                Sugar_Refer_Dict[current_seq_index]=frag_info_list[i]
        count_use_pho=0
        for k in range(len(Sugar_Only_Frag_List)):
            use_pho_model=0
            chain_id,current_seq_index,atom_name, cur_sugar_position,nuc_type,avg_score = Sugar_Only_Frag_List[k]
            if current_seq_index in Sugar_Refer_Dict:
                chain_id, current_seq_index,atom_name2, cur_pho_position,nuc_type,avg_score = Sugar_Refer_Dict[current_seq_index]
                if k+1<len(Sugar_Only_Frag_List):
                    _, next_seq_index,_, next_sugar_position,_,_ = Sugar_Only_Frag_List[k+1]
                    if next_seq_index in Sugar_Refer_Dict:
                        _, next_seq_index,_, next_pho_position,_,_ = Sugar_Refer_Dict[next_seq_index]
                        use_pho_model=1#got two consecutive phosphate, then we can use phosphate for modeling
                    else:
                        #use sugar change direction to add another virtual pho location here
                        next_sugar_position = np.array(next_sugar_position)
                        cur_sugar_position = np.array(cur_sugar_position)
                        cur_pho_position = np.array(cur_pho_position)
                        sugar_move_direction = next_sugar_position-cur_sugar_position
                        next_pho_position = sugar_move_direction+cur_pho_position
                        use_pho_model=2
                else:
                    _, prev_seq_index,_, prev_sugar_position,_,_ = Sugar_Only_Frag_List[k-1]
                    if prev_seq_index in Sugar_Refer_Dict:
                        _, prev_seq_index,_, prev_pho_position,_,_ = Sugar_Refer_Dict[prev_seq_index]
                        prev_pho_position = np.array(prev_pho_position)
                        cur_sugar_position = np.array(cur_sugar_position)
                        cur_pho_position = np.array(cur_pho_position)
                        pho_move_direction = cur_pho_position-prev_pho_position
                        next_pho_position = pho_move_direction+cur_pho_position
                        use_pho_model=3
                    else:
                        prev_sugar_position = np.array(prev_sugar_position)
                        cur_sugar_position = np.array(cur_sugar_position)
                        cur_pho_position = np.array(cur_pho_position)
                        sugar_move_direction = cur_sugar_position-prev_sugar_position
                        next_pho_position = sugar_move_direction+cur_pho_position
                        use_pho_model=4
            else:
                #if we can find previous pho/next pho, we also used pho to model
                if k>=1 and k+1<len(Sugar_Only_Frag_List):
                    _, prev_seq_index,_, prev_sugar_position,_,_ = Sugar_Only_Frag_List[k-1]
                    _, next_seq_index,_, next_sugar_position,_,_ = Sugar_Only_Frag_List[k+1]

                    if prev_seq_index in Sugar_Refer_Dict and next_seq_index in Sugar_Refer_Dict:

                        _, prev_seq_index,_, prev_pho_position,_,_ = Sugar_Refer_Dict[prev_seq_index]
                        _, next_seq_index,_, next_pho_position,_,_ = Sugar_Refer_Dict[next_seq_index]

                        prev_pho_position = np.array(prev_pho_position)
                        next_pho_position = np.array(next_pho_position)
                        cur_pho_position = (prev_pho_position+next_pho_position)/2
                        use_pho_model=5
                    elif next_seq_index in Sugar_Refer_Dict:
                        _, next_seq_index,_, next_pho_position,_,_ = Sugar_Refer_Dict[next_seq_index]
                        next_sugar_position = np.array(next_sugar_position)
                        cur_sugar_position = np.array(cur_sugar_position)
                        next_pho_position = np.array(next_pho_position)
                        sugar_move_direction = next_sugar_position-cur_sugar_position
                        cur_pho_position = next_pho_position-sugar_move_direction
                        use_pho_model=6

                    elif prev_seq_index in Sugar_Refer_Dict:
                        _, prev_seq_index,_, prev_pho_position,_,_ = Sugar_Refer_Dict[prev_seq_index]
                        next_sugar_position = np.array(next_sugar_position)
                        prev_sugar_position = np.array(prev_sugar_position)
                        cur_sugar_position = np.array(cur_sugar_position)
                        prev_pho_position = np.array(prev_pho_position)


                        sugar_move_direction = cur_sugar_position-prev_sugar_position
                        cur_pho_position = prev_pho_position+sugar_move_direction
                        sugar_move_direction = next_sugar_position-cur_sugar_position
                        next_pho_position = cur_pho_position+sugar_move_direction
                        use_pho_model=7
                    else:
                        use_pho_model=0#no way to build structure via sugar
                elif k==0:
                    _, next_seq_index,_, next_sugar_position,_,_ = Sugar_Only_Frag_List[k+1]
                    if next_seq_index in Sugar_Refer_Dict:
                        _, next_seq_index,_, next_pho_position,_,_ = Sugar_Refer_Dict[next_seq_index]
                        next_sugar_position = np.array(next_sugar_position)
                        cur_sugar_position = np.array(cur_sugar_position)
                        next_pho_position = np.array(next_pho_position)
                        sugar_move_direction = next_sugar_position-cur_sugar_position
                        cur_pho_position = next_pho_position-sugar_move_direction
                        use_pho_model=8
                    else:
                        use_pho_model=0
                elif (k+1)==len(Sugar_Only_Frag_List):
                    _, prev_seq_index,_, prev_sugar_position,_,_ = Sugar_Only_Frag_List[k-1]
                    if prev_seq_index in Sugar_Refer_Dict:
                        _, prev_seq_index,_, prev_pho_position,_,_ = Sugar_Refer_Dict[prev_seq_index]
                        prev_sugar_position = np.array(prev_sugar_position)
                        cur_sugar_position = np.array(cur_sugar_position)
                        prev_pho_position = np.array(prev_pho_position)
                        sugar_move_direction = cur_sugar_position-prev_sugar_position
                        cur_pho_position = prev_pho_position+sugar_move_direction
                        next_pho_position = cur_pho_position+sugar_move_direction
                        use_pho_model=9
                    else:
                        prev_sugar_position = np.array(prev_sugar_position)
                        cur_sugar_position = np.array(cur_sugar_position)

                        sugar_move_direction = cur_sugar_position-prev_sugar_position
                        next_sugar_position = cur_sugar_position+sugar_move_direction
                        use_pho_model=0
            cur_base_location = refer_sugar_base_location[k]
            if nuc_type=='T' and DNA_Label is False:
                nuc_type='U'
            if nuc_type=='U' and DNA_Label is True:
                nuc_type='T'
            if DNA_Label is True:
                nuc_type ="D"+nuc_type

            #try:
            if use_pho_model>0:
                print("use pho situation %d"%use_pho_model)
                count_use_pho+=1
                #assign structures based on
                atom_list,coordinate_list = build_base_atomic_nuc_rigid(cur_pho_position,next_pho_position, cur_base_location,nuc_type,map_info_list)
            else:
                atom_list,coordinate_list = build_base_atomic_nuc_rigid_sugar(cur_sugar_position,next_sugar_position, cur_base_location,nuc_type,map_info_list)
            # except:
            #     atom_list = ["C4'"]
            #     coordinate_list = [cur_sugar_position]
            #     print("svd solution failed for all atom building for this residue")
            for k in range(len(atom_list)):
                line =""
                line += "ATOM %-10d C %-3s . %-3s %-2s . %-10d  .  " % (Natm, atom_list[k], nuc_type, chain_id, current_seq_index)
                line += "%-8.3f %-8.3f %-8.3f %-6.2f %-6.8f %-10d %-2s 1 \n" % (
                coordinate_list[k][0], coordinate_list[k][1], coordinate_list[k][2], 1.0,avg_score,Natm,chain_id)
                file.write(line)
                with open(pdb_save_path,"a+") as wfile:
                    line=""
                    if len(chain_id)>2:
                        chain_id=chain_id[:2]#for saving purposes
                    line += "ATOM%7d %-4s %3s%2s%4d    " % (Natm, atom_list[k],nuc_type, chain_id,current_seq_index)
                     # tmp_dens=(point.merged_cd_dens[i,3]-dens_min)/(dens_max-dens_min)
                    line = line + "%8.3f%8.3f%8.3f%6.2f%6.2f\n" % (
                     coordinate_list[k][0], coordinate_list[k][1], coordinate_list[k][2], 1.0, avg_score)
                    wfile.write(line)

                Natm += 1
        print("we have %d/%d atomic build upon pho atoms"%(count_use_pho,len(Sugar_Only_Frag_List)))
