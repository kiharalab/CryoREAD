
from collections import defaultdict
import numpy as np
import os
from ops.os_operation import mkdir
from graph.io_utils import append_cif_info
from atomic.io_utils import Write_Atomic_Fraginfo_cif
from graph.structure_utils import clean_pho_assign_info_list
from graph.LDP_ops import Convert_LDPcoord_To_Reallocation
from ops.cif_utils import cif2pdb
def build_cif_PS(current_path_dict,path_id,ldp_sugar_location,
        Path_P_align_list,pho_point,map_info_list):
    pho_merged_cd = pho_point.merged_cd_dens[:,:3]
    current_seq_info = current_path_dict
    match_seq = current_seq_info['match_seq']
    cur_score = current_seq_info['score']
    current_ldp_location = ldp_sugar_location
    current_pho_ldp_align_info = Path_P_align_list
    current_pho_ldp_location  = [pho_merged_cd[int(kk)] for kk in current_pho_ldp_align_info ]
    current_pho_ldp_location = Convert_LDPcoord_To_Reallocation(current_pho_ldp_location , map_info_list)
    #start to ensemble all the information
    fragment_info_list = []
    current_seq_index = 1
    if path_id>=62:
        path_id= path_id%62
    chain_dict={}
    for k in range(26):
        chain_dict[k]=chr(65+k)
    for k in range(26):
        chain_dict[26+k]=chr(97+k)
    for k in range(10):
        chain_dict[52+k]="%d"%k
    align_length= len([x for x in match_seq if x!="-"])
    avg_score = cur_score/align_length
    chain_id = chain_dict[int(path_id)]
    for k in range(len(match_seq)):
        if match_seq[k]=="-":
            continue
        else:
            cur_pho_tmp_position = current_pho_ldp_location[k]
            if current_pho_ldp_align_info[k]>=0:#make sure it has assigned P
                fragment_info_list.append([chain_id,current_seq_index+1,'P',cur_pho_tmp_position,match_seq[k],avg_score])
            cur_resi_position = current_ldp_location[k]
            fragment_info_list.append([chain_id,current_seq_index+1,"C4'",cur_resi_position,match_seq[k],avg_score])

            current_seq_index+=1

    return fragment_info_list,cur_score



def build_atomic_fragment_cluster_cif_SP(Path_Assign_Dict,all_sugar_location,Path_P_align_list,Path_P_reverse_align_list,
                 save_dir,pho_point,map_info_list,refer_base_location,ext_name="geo"):
    overall_score=0
    check_file_list=[]
    for ldp_path_id in Path_Assign_Dict.keys():
        cur_ldp_sugar_location= all_sugar_location[ldp_path_id]
        cur_pho_ldp_nodeid_list = Path_P_align_list[ldp_path_id]
        cur_pho_ldp_nodeid_reverse_list = Path_P_reverse_align_list[ldp_path_id]
        current_path_dict = Path_Assign_Dict[ldp_path_id]
        if current_path_dict['direction']==1:
            input_pho_nodeid = cur_pho_ldp_nodeid_list
        else:
            input_pho_nodeid = cur_pho_ldp_nodeid_reverse_list
        tmp_save_path = os.path.join(save_dir,"%s_%d_path.cif"%(ext_name, ldp_path_id))

        fragment_info_list,frag_score = build_cif_PS(current_path_dict,ldp_path_id, cur_ldp_sugar_location,
            input_pho_nodeid,pho_point,map_info_list)
        overall_score+=frag_score
        if len(fragment_info_list)<=2:
            #not enough to build an atomic model, skip it.
            continue
        #clean pho locations in case we have one pho assigned to 2 sugars in the path
        fragment_info_list,_ = clean_pho_assign_info_list(fragment_info_list,pho_point.merged_cd_dens[:,:3],map_info_list)
        fragment_info_list,further_flag = clean_pho_assign_info_list(fragment_info_list,pho_point.merged_cd_dens[:,:3],map_info_list,round=2)
        if further_flag:
            fragment_info_list,_ = clean_pho_assign_info_list(fragment_info_list,pho_point.merged_cd_dens[:,:3],map_info_list,round=2)
        if len(fragment_info_list)<=2:
            #not enough to build an atomic model, skip it.
            continue
        Write_Atomic_Fraginfo_cif("%s_%d_cluster"%(ext_name,ldp_path_id),fragment_info_list,refer_base_location,tmp_save_path,False,map_info_list)
        check_file_list.append(tmp_save_path)

    return check_file_list,overall_score

def Build_Atomic_Structure(overall_dict,
    all_sugar_location, save_dir,
    Path_P_align_list,Path_P_reverse_align_list,pho_point,map_info_list,refer_base_location):
    frag_dir = os.path.join(save_dir,"Path_Atomic_Frags")
    mkdir(frag_dir)

    check_file_list,overall_score = build_atomic_fragment_cluster_cif_SP(overall_dict,
        all_sugar_location,Path_P_align_list,Path_P_reverse_align_list,
                 frag_dir,pho_point,map_info_list,refer_base_location)

    fragment_all_path = os.path.join(save_dir,"Final_Assemble_geo.cif")
    with open(fragment_all_path,'w') as file:
        file.write("#score: %f\n"%overall_score)
    check_file_list.sort()
    for item in check_file_list:
        #cur_frag_path = os.path.join(save_dir,item)
        cur_frag_name = os.path.split(item)[1]
        cur_entry_id = cur_frag_name.replace(".cif","")
        append_cif_info(cur_entry_id,item,fragment_all_path)
    fragment_all_path = os.path.join(save_dir,"Final_Assemble_geo.pdb")
    with open(fragment_all_path,'w') as file:
        file.write("#score: %f\n"%overall_score)
    check_file_list.sort()
    pdb_file_list=[]
    for item in check_file_list:
        pdb_file_name = item.replace(".cif",".pdb")
        cif2pdb(item,pdb_file_name)
        pdb_file_list.append(pdb_file_name)
    with open(fragment_all_path,'a+') as wfile:
        for item in pdb_file_list:
            with open(item,'r') as rfile:
                for line in rfile:
                    wfile.write(line)
        

