
from collections import defaultdict
import numpy as np
import os
from ops.os_operation import mkdir
from graph.io_utils import append_cif_info
from ops.cif_utils import cif2pdb
from graph.structure_utils import build_noclusters_extra,build_clusters,\
    build_clusters_advanced,merge_cluster_cif_PS,clean_pho_assign_info_list
from atomic.io_utils import Write_Atomic_Fraginfo_cif


def build_atomic_fragment_cluster_cif_SP(Path_Assign_Dict,all_sugar_location,Path_P_align_list,Path_P_reverse_align_list,
                 ldp_size,save_dir,pho_point,map_info_list,refer_base_location,DNA_Label,ext_name="path"):
    overall_score=0
    check_file_list=[]

    for ldp_path_id in Path_Assign_Dict.keys():
        cur_ldp_sugar_location= all_sugar_location[ldp_path_id]
        cur_pho_ldp_nodeid_list = Path_P_align_list[ldp_path_id]
        cur_pho_ldp_nodeid_reverse_list = Path_P_reverse_align_list[ldp_path_id]
        current_path_dict = Path_Assign_Dict[ldp_path_id]
        all_starting_index_keys = [int(x) for x in current_path_dict.keys()]#should already be in order
        all_seq_assign_length = [len(current_path_dict[tmp_path_id]['match_seq']) for tmp_path_id in current_path_dict]
        assert len(all_starting_index_keys)==len(np.unique(all_starting_index_keys))
        all_starting_index_keys = np.array(all_starting_index_keys)
        all_seq_assign_length = np.array(all_seq_assign_length)
        sorted_indexes = np.argsort(all_starting_index_keys)
        all_starting_index_keys = all_starting_index_keys[sorted_indexes]
        all_seq_assign_length = all_seq_assign_length[sorted_indexes]
        if 'extra' in ext_name:
            interval_clusters = build_noclusters_extra(all_starting_index_keys,all_seq_assign_length)
        else:
            interval_clusters = build_clusters_advanced(all_starting_index_keys,all_seq_assign_length)#overlaped regions merged as one framenet. non-overlap will be different fragments
        print("we clustered %d clusters out of %d fragments"%(len(interval_clusters),len(all_starting_index_keys)))
        for cluster_id in interval_clusters:
            current_cluster_id_list = interval_clusters[cluster_id]
            tmp_save_path = os.path.join(save_dir,"%s_%d_cluster_%d.cif"%(ext_name, ldp_path_id,cluster_id))
            tmp_verify_path = os.path.join(save_dir,"%s_%d_cluster_%d_verify.txt"%(ext_name,ldp_path_id,cluster_id))
            fragment_info_list,frag_score = merge_cluster_cif_PS(current_path_dict,current_cluster_id_list,cur_ldp_sugar_location,
            tmp_verify_path,cur_pho_ldp_nodeid_list,cur_pho_ldp_nodeid_reverse_list,pho_point,map_info_list)
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
            Write_Atomic_Fraginfo_cif("%s_%d_cluster_%d"%(ext_name,ldp_path_id,cluster_id),fragment_info_list,refer_base_location,tmp_save_path,DNA_Label,map_info_list)
            check_file_list.append(tmp_save_path)

    return check_file_list,overall_score

def Build_Atomic_Structure(solve_frag_combine_list,order_key_index,order_chain_index,overall_dict,
    all_sugar_location, ldp_size,save_dir,checking_stride,top_select,
    Path_P_align_list,Path_P_reverse_align_list,pho_point,map_info_list,Extra_Added_Assign_Dict,refer_base_location,DNA_Label):
    Path_Assign_Dict=defaultdict(dict)#[path_id][assignment starting_index]:[seq_info]
    for order_index in solve_frag_combine_list:
        order_index = int(order_index)
        current_key1 = order_key_index[order_index]
        current_chain_candidate1 = order_chain_index[order_index]
        split_key1 = current_key1.split("_")
        ldp_path1= int(split_key1[0])
        ldp_starting_index1 = int(split_key1[1])
        current_seq_info1 = overall_dict[current_key1][current_chain_candidate1]
        Path_Assign_Dict[ldp_path1][ldp_starting_index1]=current_seq_info1
    #save all non-overlap in one dir, while extra in another dir.
    non_overlap_dir = os.path.join(save_dir,"CP_SAT_frags")
    mkdir(non_overlap_dir)
    check_file_list,overall_score1 = build_atomic_fragment_cluster_cif_SP(Path_Assign_Dict,all_sugar_location,Path_P_align_list,Path_P_reverse_align_list,
                 ldp_size,non_overlap_dir,pho_point,map_info_list,refer_base_location,DNA_Label,ext_name="path")


    overlap_dir = os.path.join(save_dir,"extra_support_frags")
    mkdir(overlap_dir)
    check_file_list2,overall_score2 = build_atomic_fragment_cluster_cif_SP(Extra_Added_Assign_Dict,all_sugar_location,Path_P_align_list,Path_P_reverse_align_list,
                 ldp_size,overlap_dir,pho_point,map_info_list,refer_base_location,DNA_Label,ext_name="extra")

    check_file_list = check_file_list+check_file_list2

    overall_score = overall_score1+overall_score2
    fragment_all_path = os.path.join(save_dir,"Final_Assemble_%d_%d_%d.cif"%(ldp_size,checking_stride,top_select))
    with open(fragment_all_path,'w') as file:
        file.write("#score: %f\n"%overall_score)
    check_file_list.sort()
    for item in check_file_list:
        #cur_frag_path = os.path.join(save_dir,item)
        cur_frag_name = os.path.split(item)[1]
        cur_entry_id = cur_frag_name.replace(".cif","")
        append_cif_info(cur_entry_id,item,fragment_all_path)

    fragment_all_path = os.path.join(save_dir,"Final_Assemble_%d_%d_%d.pdb"%(ldp_size,checking_stride,top_select))
    with open(fragment_all_path,'w') as file:
        file.write("#score: %f\n"%overall_score)
    check_file_list.sort()
    Natm=1
    Nres=0
    for item in check_file_list:
        cur_pdb_name = item[:-4]+".pdb"
        prev_resi = None
        with open(cur_pdb_name,'r') as rfile:
            with open(fragment_all_path,"a+") as wfile:
                for line in rfile:
                    #advance modifying the residue id to build

                    if (line.startswith('ATOM')):
                        #chain_id, current_seq_index,atom_name2, cur_pho_position,nuc_type,avg_score
                        chain_name = line[21]
                        atom_name = line[12:16]
                        x=float(line[30:38])
                        y=float(line[38:46])
                        z=float(line[46:55])
                        resi=int(line[22:26])
                        score = float(line[60:68])
                        resn = line[17:20]
                        if resi!=prev_resi:
                            Nres+=1
                            prev_resi=resi
                        line=""
                        line += "ATOM%7d %-4s %3s%2s%4d    " % (Natm, atom_name,resn, chain_name,Nres)
                        line = line + "%8.3f%8.3f%8.3f%6.2f%6.2f\n" % (x,y,z, 1.0, score)
                        wfile.write(line)

                        Natm+=1
                    #wfile.write(line)
                    else:
                        wfile.write(line)



def Build_Atomic_Model_nonoverlap_frag(Path_Assign_Dict,
    all_sugar_location, ldp_size,save_dir,checking_stride,top_select,
    Path_P_align_list,Path_P_reverse_align_list,pho_point,map_info_list,
    Extra_Added_Assign_Dict,refer_base_location,DNA_Label):
    non_overlap_dir = os.path.join(save_dir,"CP_SAT_frags")
    mkdir(non_overlap_dir)
    check_file_list,overall_score1 = build_atomic_fragment_cluster_cif_SP(Path_Assign_Dict,all_sugar_location,Path_P_align_list,Path_P_reverse_align_list,
                 ldp_size,non_overlap_dir,pho_point,map_info_list,refer_base_location,DNA_Label,ext_name="extra")


    overlap_dir = os.path.join(save_dir,"extra_support_frags")
    mkdir(overlap_dir)
    check_file_list2,overall_score2 = build_atomic_fragment_cluster_cif_SP(Extra_Added_Assign_Dict,all_sugar_location,Path_P_align_list,Path_P_reverse_align_list,
                 ldp_size,overlap_dir,pho_point,map_info_list,refer_base_location,DNA_Label,ext_name="extra2")

    check_file_list = check_file_list+check_file_list2

    overall_score = overall_score1+overall_score2
    fragment_all_path = os.path.join(save_dir,"Final_Assemble_%d_%d_%d.cif"%(ldp_size,checking_stride,top_select))
    with open(fragment_all_path,'w') as file:
        file.write("#score: %f\n"%overall_score)
    check_file_list.sort()
    for item in check_file_list:
        #cur_frag_path = os.path.join(save_dir,item)
        cur_frag_name = os.path.split(item)[1]
        cur_entry_id = cur_frag_name.replace(".cif","")
        append_cif_info(cur_entry_id,item,fragment_all_path)

    fragment_all_path = os.path.join(save_dir,"Final_Assemble_%d_%d_%d.pdb"%(ldp_size,checking_stride,top_select))
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
