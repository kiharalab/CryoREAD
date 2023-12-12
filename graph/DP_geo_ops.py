


from collections import defaultdict
from turtle import screensize
import numpy as np
from data_processing.map_utils import permute_ns_coord_to_pdb,process_map_data,permute_pdb_coord_to_map,permute_map_coord_to_pdb
from ops.math_calcuation import calculate_distance
import os
from ops.os_operation import mkdir
import pickle
from graph.LDP_ops import Convert_LDPcoord_To_Reallocation
from numba import jit

def Calculate_Distance_array(fragment_ldp_location):
    distance_array =np.zeros([len(fragment_ldp_location),len(fragment_ldp_location)])
    for i in range(len(fragment_ldp_location)):
        for j in range(i+1,len(fragment_ldp_location)):
            distance_array[i,j] = calculate_distance(fragment_ldp_location[i],fragment_ldp_location[j])
            distance_array[j,i]= distance_array[i,j]
    return distance_array
@jit(nogil=True,nopython=True)
def global_align_score(ldp_gap_penalty,match_matrix,
    pointer, scratch,ldp_distance_array, ldp_prev_connect):
    low_cut_off_distance=4
    high_cut_off_distance=10
    mean_distance=7
    #two operations for pointer: 1 use, 0 skip
    for i in range(len(ldp_distance_array)):
        if i==0:
            pointer[i]=1
            scratch[i]=match_matrix[i]
            ldp_prev_connect[i]=-1 #indicates itself is the last connected
        else:
            #1 choose
            prev_neighbor = int(ldp_prev_connect[i-1]) if ldp_prev_connect[i-1]!=-1 else i-1
            #choose_best = scratch[prev_neighbor]+match_matrix[i]
            cur_distance = ldp_distance_array[prev_neighbor,i]

            if cur_distance<=low_cut_off_distance or cur_distance>=high_cut_off_distance:
                current_best= scratch[i-1]+match_matrix[i]-ldp_gap_penalty*abs(cur_distance-mean_distance)**2#ldp_gap
            else:
                current_best= scratch[i-1]+match_matrix[i]

            #2 skip this node
            reasonable_skip_ldp=0
            if cur_distance+ldp_distance_array[i,i+1]<high_cut_off_distance:
                reasonable_skip_ldp=1
            if  cur_distance<=low_cut_off_distance:
                reasonable_skip_ldp=1
            if reasonable_skip_ldp:
                left_cur = scratch[i-1]
            else:
                left_cur = scratch[i-1]-ldp_gap_penalty*abs(cur_distance-mean_distance)**2#ldp_gap
            if current_best>=left_cur:
                scratch[i]=current_best
                pointer[i]=1
                ldp_prev_connect[i]=-1
            else:
                scratch[i]= left_cur
                pointer[i]=0
                ldp_prev_connect[i]=prev_neighbor
    return scratch,pointer,ldp_prev_connect

def dynamic_assign_geo(updated_base_list,ldp_gap_penalty,fragment_distance_array,save_path=None):
    map_dict={0:"A",1:"T",2:"C",3:"G"}
    match_matrix = np.max(updated_base_list,axis=1)#updated_base_list
    if save_path is not None:
        score_path = os.path.join(save_path,"match_score.txt")
        np.savetxt(score_path,match_matrix)
    pointer = np.zeros(len(updated_base_list))#indicates the operations at position i
    scratch = np.zeros(len(updated_base_list))
    ldp_prev_connect = np.zeros(len(updated_base_list))-1
    scratch,pointer,ldp_prev_connect=global_align_score(ldp_gap_penalty,match_matrix,
    pointer, scratch,fragment_distance_array, ldp_prev_connect)
    if save_path is not None:
        score_path = os.path.join(save_path,"optimal_score.txt")
        np.savetxt(score_path,scratch)
        score_path = os.path.join(save_path,"optimal_direction.txt")
        np.savetxt(score_path,pointer)#verfied corret
        score_path = os.path.join(save_path,"ldp_prev_connect.txt")
        np.savetxt(score_path,ldp_prev_connect)
    input_seq_line=""
    for i in range(len(updated_base_list)):
        choice=int(np.argmax(updated_base_list[i]))
        input_seq_line += map_dict[choice]
    # if len(scratch)>5:
    #     max_score_index = np.argmax(scratch[-5:])+len(scratch)-5
    # else:
    #     max_score_index =  np.argmax(scratch)
    # match_matrix = np.zeros(len(updated_base_list))-1
    # max_score_index = int(max_score_index)
    # match_matrix[max_score_index]=1
    # prev_index = pointer[int(max_score_index)]
    # while check_pointer!=0:
    #for any ldp_prev_connect=-1, we select it.
    max_score = np.max(scratch)
    match_matrix = np.zeros(len(updated_base_list))
    match_matrix[ldp_prev_connect==-1]=1
    if save_path is not None:
        score_path = os.path.join(save_path,"match_seq.txt")
        match_seq=""
        with open(score_path,'w') as file:
            file.write(input_seq_line+"\n")
            for k in range(len(match_matrix)):
                if match_matrix[k]==1:
                    file.write(input_seq_line[k])
                    match_seq +=input_seq_line[k]
                else:
                    file.write("-")
                    match_seq+="-"
            file.write("\n")
    else:
        match_seq=""
        for k in range(len(match_matrix)):
            if match_matrix[k]==1:
                match_seq +=input_seq_line[k]
            else:
                match_seq+="-"
    return max_score,match_seq







def greedy_assign_geo(All_Base_Path_List_sugar,All_Path_List_sugar,Path_P_align_list,Path_P_reverse_align_list,
        Pho_Prob_Refer_Dict,Pho_Prob_Refer_Reverse_Dict,
     sugar_point,map_info_list,greedy_save_path):
    cur_frag_path = os.path.join(greedy_save_path,"All_fragment_geo.pkl")
    cur_frag_path2 = os.path.join(greedy_save_path,"Location_fragment_geo.pkl")
    if os.path.exists(cur_frag_path) and os.path.exists(cur_frag_path2):
        with open(cur_frag_path, 'rb') as handle:
            overall_dict= pickle.load(handle)
        with open(cur_frag_path2, 'rb') as handle:
            frag_location_dict = pickle.load(handle)
        return overall_dict,frag_location_dict
    merged_cd_dens = sugar_point.merged_cd_dens
    ldp_proper_range = [4,10]
    overall_dict={}
    frag_location_dict = {}

    ldp_gap_penalty = 10 # only applied if not reasonable skip
    for path_id,cur_path_list in enumerate(All_Base_Path_List_sugar):
        cur_path_save_path = os.path.join(greedy_save_path,"path_%d"%path_id)
        mkdir(cur_path_save_path)
        current_path = All_Path_List_sugar[path_id]
        current_pho_align = Path_P_align_list[path_id]
        current_pho_align_reverse = Path_P_reverse_align_list[path_id]
        current_length = len(cur_path_list)
        current_base_list = All_Base_Path_List_sugar[path_id]
        updated_base_list = [np.array(current_base_list[j]) for j in range(len(current_base_list))]
        updated_base_list = np.array(updated_base_list)
        #update the reward matrix based on the pho updates
        pho_order_base_prob_list =[np.array(Pho_Prob_Refer_Dict[int(k)]) if k!=-1 and int(k) in Pho_Prob_Refer_Dict  else np.zeros(4) for k in current_pho_align]
        pho_reverse_base_prob_list = [np.array(Pho_Prob_Refer_Reverse_Dict[int(k)]) if k!=-1 and int(k) in Pho_Prob_Refer_Reverse_Dict else np.zeros(4) for k in current_pho_align_reverse]
        assert len(pho_order_base_prob_list)==len(updated_base_list) and len(pho_order_base_prob_list)==len(pho_reverse_base_prob_list)
        updated_base_score_list = [updated_base_list[j]+pho_order_base_prob_list[j]+pho_reverse_base_prob_list[j] for j in range(len(current_base_list))]
        updated_base_score_list = np.array(updated_base_score_list)
        study_ldp_base_list = updated_base_score_list

        #study_ldp_base_list_reverse = updated_base_score_reverse_list[start_ldp_index:end_ldp_index]
        current_location_list = [merged_cd_dens[int(kk)] for kk in current_path]
        dp_save_path = os.path.join(cur_path_save_path,"path_order")
        mkdir(dp_save_path)
        fragment_ldp_location= Convert_LDPcoord_To_Reallocation(current_location_list, map_info_list)
        fragment_distance_array = Calculate_Distance_array(fragment_ldp_location)
        max_score, match_seq_line = dynamic_assign_geo(study_ldp_base_list,ldp_gap_penalty,
        fragment_distance_array, dp_save_path)
        dp_save_path = os.path.join(cur_path_save_path,"path_reverse_order")
        mkdir(dp_save_path)
        max_score2, match_seq_line2 = dynamic_assign_geo(study_ldp_base_list[::-1],ldp_gap_penalty,
        fragment_distance_array[::-1,::-1], dp_save_path)
        new_dict={}
        if max_score>max_score2:
            match_seq_line=match_seq_line
            new_dict['direction']=1
        else:
            match_seq_line= match_seq_line2[::-1]
            new_dict['direction']=-1
        new_dict['score']=max_score
        new_dict['match_seq']=match_seq_line
        overall_dict[path_id]=new_dict
        frag_location_dict[path_id]=fragment_ldp_location
    print("geometrical based dp finished selecting ldp for atomic modeling")
    with open(cur_frag_path, 'wb') as handle:
        pickle.dump(overall_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    with open(cur_frag_path2, 'wb') as handle:
        pickle.dump(frag_location_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return overall_dict,frag_location_dict

