

from collections import defaultdict
import numpy as np
import os
from ops.os_operation import mkdir
from graph.DP_ops import Calculate_Distance_array,dynamic_assign_multi
from graph.structure_modeling import build_clusters_advanced
import pickle


def reassign_basedonfrag(solve_frag_combine_list,order_key_index,order_chain_index,overall_dict,
    all_sugar_location, ldp_size,save_dir,checking_stride,top_select,chain_dict,
    All_Base_Path_List_sugar,All_Path_List_sugar,Path_P_align_list,Path_P_reverse_align_list,
        Pho_Prob_Refer_Dict,Pho_Prob_Refer_Reverse_Dict,):
    cur_frag_path = os.path.join(save_dir,"Fragment_reassign_%d_%d_%d.pkl"%(ldp_size,checking_stride,top_select))
    if os.path.exists(cur_frag_path):
        with open(cur_frag_path, 'rb') as handle:
            Path_Assign_Dict= pickle.load(handle)

        return Path_Assign_Dict
    Path_Assign_Dict=defaultdict(dict)#[path_id][assignment starting_index]:[seq_info]
    ldp_gap_penalty = 10 # only applied if not reasonable skip
    seq_gap_penalty = 25 # applied once we need to skip an amino acid
    for order_index in solve_frag_combine_list:
        order_index = int(order_index)
        current_key1 = order_key_index[order_index]
        current_chain_candidate1 = order_chain_index[order_index]
        split_key1 = current_key1.split("_")
        ldp_path1= int(split_key1[0])
        ldp_starting_index1 = int(split_key1[1])
        current_seq_info1 = overall_dict[current_key1][current_chain_candidate1]
        Path_Assign_Dict[ldp_path1][ldp_starting_index1]=current_seq_info1
    old_frag_path = os.path.join(save_dir,"Fragment_assemble_%d_%d_%d.pkl"%(ldp_size,checking_stride,top_select))
    with open(old_frag_path, 'wb') as handle:
        pickle.dump(Path_Assign_Dict, handle, protocol=pickle.HIGHEST_PROTOCOL)


    Path_ReAssign_Dict=defaultdict(dict)
    for ldp_path_id in Path_Assign_Dict.keys():
        cur_ldp_sugar_location= all_sugar_location[ldp_path_id]
        current_path_dict = Path_Assign_Dict[ldp_path_id]
        all_starting_index_keys = [int(x) for x in current_path_dict.keys()]#should already be in order
        all_seq_assign_length = [len(current_path_dict[tmp_path_id]['match_seq']) for tmp_path_id in current_path_dict]
        assert len(all_starting_index_keys)==len(np.unique(all_starting_index_keys))
        all_starting_index_keys = np.array(all_starting_index_keys)
        all_seq_assign_length = np.array(all_seq_assign_length)
        sorted_indexes = np.argsort(all_starting_index_keys)
        all_starting_index_keys = all_starting_index_keys[sorted_indexes]
        all_seq_assign_length = all_seq_assign_length[sorted_indexes]
        interval_clusters = build_clusters_advanced(all_starting_index_keys,all_seq_assign_length)
        current_path = All_Path_List_sugar[ldp_path_id]
        current_pho_align = Path_P_align_list[ldp_path_id]
        current_pho_align_reverse = Path_P_reverse_align_list[ldp_path_id]

        current_base_list = All_Base_Path_List_sugar[ldp_path_id]
        updated_base_list = [np.array(current_base_list[j]) for j in range(len(current_base_list))]
        updated_base_list = np.array(updated_base_list)
        #update the reward matrix based on the pho updates
        pho_order_base_prob_list =[np.array(Pho_Prob_Refer_Dict[int(k)]) if k!=-1 and int(k) in Pho_Prob_Refer_Dict  else np.zeros(4) for k in current_pho_align]
        pho_reverse_base_prob_list = [np.array(Pho_Prob_Refer_Reverse_Dict[int(k)]) if k!=-1 and int(k) in Pho_Prob_Refer_Reverse_Dict else np.zeros(4) for k in current_pho_align_reverse]

        assert len(pho_order_base_prob_list)==len(updated_base_list) and len(pho_order_base_prob_list)==len(pho_reverse_base_prob_list)
        updated_base_score_list = [updated_base_list[j]+pho_order_base_prob_list[j]+pho_reverse_base_prob_list[j] for j in range(len(current_base_list))]
        updated_base_score_list = np.array(updated_base_score_list)
        for cluster_id in interval_clusters:
            current_cluster_id_list = interval_clusters[cluster_id]
            current_cluster_id_list.sort()

            begin_index = current_cluster_id_list[0]
            if len(current_cluster_id_list)==1:
                #for non-overlapped one, simply take them.
                Path_ReAssign_Dict[ldp_path_id][begin_index]=current_path_dict[begin_index]
                continue
            #cur_ldp_size = len(current_path_dict[current_cluster_id_list[-1]]['match_seq'])
            #end_index = current_cluster_id_list[-1]+cur_ldp_size
            #check all the matches to get the final extended region
            end_index=-1
            for k in range(len(current_cluster_id_list)):
                cur_ldp_size = len(current_path_dict[current_cluster_id_list[k]]['match_seq'])
                cur_end_index = current_cluster_id_list[k]+cur_ldp_size
                if end_index<=cur_end_index:
                    end_index=cur_end_index

            start_ldp_index=begin_index
            end_ldp_index=end_index
            fragment_ldp_location = cur_ldp_sugar_location[start_ldp_index:end_ldp_index]
            #begin index, end_index suggested the ldps that needs for current re-dp.

            #next step: find the seqeucen to assign.
            current_seq_info = current_path_dict[current_cluster_id_list[0]]
            overall_direction = current_seq_info['direction']
            overall_chain = current_seq_info['chain']
            overall_chain_length = current_seq_info['chain_length']
            init_interval = current_path_dict[current_cluster_id_list[0]]['interval']
            for path_id in current_cluster_id_list:
                current_interval = current_path_dict[path_id]['interval']
                if current_interval[0]<=init_interval[0]:
                    init_interval[0]=current_interval[0]
                if current_interval[1]>=init_interval[1]:
                    init_interval[1]=current_interval[1]
            print("current interval:",init_interval)

            current_chain_sequence = chain_dict[overall_chain]
            if overall_direction==-1:
                current_chain_sequence = current_chain_sequence[::-1]
            begin_seq_index = max(0,int(init_interval[0]-ldp_size/2))
            end_seq_index = min(overall_chain_length,int(init_interval[1]+ldp_size/2))
            current_chain_sequence = current_chain_sequence[begin_seq_index:end_seq_index]




            #then the fragment information to calculate dp

            study_ldp_base_list = updated_base_score_list[start_ldp_index:end_ldp_index]
            dp_save_path = os.path.join(save_dir,"path%d_starting%d_chain_%s"%(ldp_path_id,start_ldp_index,overall_chain))
            mkdir(dp_save_path)
            fragment_distance_array = Calculate_Distance_array(fragment_ldp_location)

            max_score, match_seq_line,match_seq_interval = dynamic_assign_multi(study_ldp_base_list,current_chain_sequence,
                        ldp_gap_penalty,seq_gap_penalty,fragment_distance_array, dp_save_path)
            path_assign_score_list =max_score
            path_assign_score_list=np.array(path_assign_score_list)
            path_score_index = np.argsort(path_assign_score_list)
            path_score_index = path_score_index[::-1]#from bigger to smaller
            for select_index in path_score_index[:1]:

                current_score = max_score[select_index]
                current_match_seq = match_seq_line[select_index]
                current_match_seq_interval = match_seq_interval[select_index]


            new_dict={}
            new_dict['match_seq']=current_match_seq
            new_dict['score']=current_score
            new_dict['interval']= [begin_seq_index+current_match_seq_interval[0],begin_seq_index+current_match_seq_interval[1]]
            new_dict['chain']=overall_chain
            new_dict['direction']=overall_direction
            new_dict['chain_length'] = overall_chain_length
            Path_ReAssign_Dict[ldp_path_id][start_ldp_index]= new_dict
    with open(cur_frag_path, 'wb') as handle:
        pickle.dump(Path_ReAssign_Dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return Path_ReAssign_Dict


def merge_assign_geo_seq(overall_geo_dict,overall_seq_dict):
    Final_Dict=defaultdict(dict)
    for ldp_path_id in overall_geo_dict.keys():
        current_dict = overall_geo_dict[ldp_path_id]
        current_direction = current_dict['direction']
        current_match_seq = current_dict['match_seq']
        #read sequence assignment information
        current_seqpath_dict = overall_seq_dict[ldp_path_id]
        all_starting_index_keys = [int(x) for x in current_seqpath_dict.keys()]#should already be in order
        all_seq_assign_length = [len(current_seqpath_dict[tmp_path_id]['match_seq']) for tmp_path_id in current_seqpath_dict]
        assert len(all_starting_index_keys)==len(np.unique(all_starting_index_keys))
        all_starting_index_keys = np.array(all_starting_index_keys)
        all_seq_assign_length = np.array(all_seq_assign_length)
        sorted_indexes = np.argsort(all_starting_index_keys)
        all_starting_index_keys = all_starting_index_keys[sorted_indexes]
        all_seq_assign_length = all_seq_assign_length[sorted_indexes]


        for starting_index in all_starting_index_keys:
            current_seq_info = current_seqpath_dict[int(starting_index)]
            seqassign_match_seq = current_seq_info['match_seq']
            #no matter direction 0 or 1, we always assign from left, then we also do not need to switch back the pho/sugar location again
            #cur_direction = current_seq_info['direction']
            prev_seq = current_match_seq
            current_match_seq = current_match_seq[:starting_index]+seqassign_match_seq+current_match_seq[starting_index+len(seqassign_match_seq):]
            print("reassign %d index of %s to %s"%(starting_index,prev_seq,current_match_seq))

        print("final seq+geo sequence %s"%current_match_seq)
        current_dict['match_seq']=current_match_seq
        Final_Dict[ldp_path_id]=current_dict
    return Final_Dict
