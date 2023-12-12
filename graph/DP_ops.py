from collections import defaultdict
import numpy as np
from data_processing.map_utils import permute_ns_coord_to_pdb,process_map_data,permute_pdb_coord_to_map,permute_map_coord_to_pdb
from ops.math_calcuation import calculate_distance
import os
from ops.os_operation import mkdir
import pickle
from graph.LDP_ops import  Convert_LDPcoord_To_Reallocation
from graph.visualize_utils import Visualize_assign_DPbase

def Calculate_Distance_array(fragment_ldp_location):
    distance_array =np.zeros([len(fragment_ldp_location),len(fragment_ldp_location)])
    for i in range(len(fragment_ldp_location)):
        for j in range(i+1,len(fragment_ldp_location)):
            distance_array[i,j] = calculate_distance(fragment_ldp_location[i],fragment_ldp_location[j])
            distance_array[j,i]= distance_array[i,j]
    return distance_array
from numba import jit
@jit(nogil=True,nopython=True)
def global_align_score(ldp_gap,seq_gap, matrix, ldp_list,
seq_list, pointer, scratch,seq_gap_count,ldp_gap_count,
ldp_distance_array,ldp_prev_connect,seq_gap_total_limit=2):
    low_cut_off_distance=4
    high_cut_off_distance=10
    mean_distance=7
    for i in range(1,len(ldp_list)+1):
        for j in range(1,len(seq_list)+1):
            reasonable_skip_ldp=0
            if i==1:
                current_best= scratch[i-1,j-1]+matrix[i,j]
            else:
                prev_neighbor = int(ldp_prev_connect[i-1,j-1]) if ldp_prev_connect[i-1,j-1]!=-1 else i-1
                cur_distance = ldp_distance_array[prev_neighbor,i]

                if cur_distance<=low_cut_off_distance or cur_distance>=high_cut_off_distance:
                    #current_best=-999999
                #elif cur_distance>=high_cut_off_distance:
                    current_best= scratch[i-1,j-1]+matrix[i,j]-ldp_gap*abs(cur_distance-mean_distance)**2#ldp_gap
                else:
                    current_best= scratch[i-1,j-1]+matrix[i,j]

            #left_best= 0
            #left_best_step=0
            if i==1:
               left_cur = scratch[i-1,j]
            else:
                prev_neighbor = int(ldp_prev_connect[i-1,j]) if ldp_prev_connect[i-1,j]!=-1 else i-1
                cur_distance = ldp_distance_array[prev_neighbor,i]
                if i!=len(seq_list) and cur_distance+ldp_distance_array[i,i+1]<high_cut_off_distance:
                    reasonable_skip_ldp=1
                if  cur_distance<=low_cut_off_distance:
                    reasonable_skip_ldp=1
                #if it's too far, do not encourage skip the ldp
                # if cur_distance>=high_cut_off_distance:
                # #     #left_cur = -999999 #do not allow skip this ldp point any more #scratch[i-1,j]-abs(cur_distance-mean_distance)*ldp_gap#-ldp_gap#-ldp_gap*max(ldp_gap_count[i-1,j]-1,0)
                #      left_cur = scratch[i-1,j]-100*abs(cur_distance-mean_distance)#ldp_gap
                # # else:
                # #elif ldp_gap_count[i-1,j]>=2:
                # #    left_cur = -999999
                # else:
                #     left_cur = scratch[i-1,j]
                    #else:
                    #    left_cur = scratch[i-1,j]-ldp_gap*(ldp_gap_count[i-1,j]+1)
                if reasonable_skip_ldp:
                    left_cur = scratch[i-1,j]
                else:
                    left_cur = scratch[i-1,j]-ldp_gap*abs(cur_distance-mean_distance)**2#ldp_gap
            left_best_step=1
            # up_best =0
            # up_best_step=0
            if seq_gap_count[i,j-1]<=seq_gap_total_limit:
                up_cur = scratch[i,j-1]-seq_gap#-seq_gap*(seq_gap_count[i,j-1]+1)#include the penalty for adding this
            else:
                up_cur = scratch[i,j-1]-seq_gap*(seq_gap_count[i,j-1]+1)#keep all possible combinations
            #up_cur = scratch[i,j-1]-seq_gap
            #up_cur=0#do not allow any skip for the sequences
            up_best_step=1
            #must take some operations
            final_best = max(current_best,left_cur,up_cur)#max(0,current_best,left_cur,up_cur)
            scratch[i,j]=final_best
            if final_best==up_cur:
                pointer[i,j]=up_best_step
                seq_gap_count[i,j] = 1+seq_gap_count[i,j-1]
                ldp_prev_connect[i,j] = ldp_prev_connect[i,j-1]
                ldp_gap_count[i,j]= ldp_gap_count[i,j-1]

            elif final_best==left_cur:
                pointer[i,j]=-left_best_step
                seq_gap_count[i,j] = seq_gap_count[i-1,j]
                if ldp_prev_connect[i-1,j]!=-1:
                    ldp_prev_connect[i,j] = ldp_prev_connect[i-1,j]
                else:
                    ldp_prev_connect[i,j] = i-1
                if reasonable_skip_ldp:
                    ldp_gap_count[i,j]=ldp_gap_count[i-1,j]
                else:
                    ldp_gap_count[i,j]=ldp_gap_count[i-1,j]+1

            else:
                pointer[i,j]=999999#indicates directly go diagnol
                seq_gap_count[i,j] = seq_gap_count[i-1,j-1]
                ldp_gap_count[i,j]= ldp_gap_count[i-1,j-1]

    return scratch,pointer,seq_gap_count,ldp_prev_connect,ldp_gap_count


def dynamic_assign_multi(updated_base_list,current_chain,
        ldp_gap_penalty,seq_gap_penalty,fragment_distance_array,save_path=None):
    map_dict={0:"A",1:"T",2:"C",3:"G"}
    #build M*N matrix
    match_matrix = np.zeros([len(updated_base_list)+1,len(current_chain)+1])
    for i in range(len(updated_base_list)):
        current_prob = updated_base_list[i]
        for j in range(len(current_chain)):
            cur_label=int(current_chain[j])
            match_matrix[i+1,j+1]=current_prob[cur_label]*100#-100#if it's higher than reference value, it should be bigger than 1
    if save_path is not None:
        score_path = os.path.join(save_path,"match_score.txt")
        np.savetxt(score_path,match_matrix)
    match_score = match_matrix
    pointer = np.zeros([len(updated_base_list)+1,len(current_chain)+1])
    scratch = np.zeros([len(updated_base_list)+1,len(current_chain)+1])
    #init the scracth matrix
    for k in range(1,len(updated_base_list)+1):
       scratch[k,0]=-999999#do not allow open gap
    #for k in range(1,len(current_chain)+1):
    #   scratch[0,k]=-seq_gap_penalty*k

    fragment_final_distance = np.zeros([len(fragment_distance_array)+1,len(fragment_distance_array)+1])
    fragment_final_distance[1:,1:]=fragment_distance_array
    fragment_final_distance[0,:]=20
    fragment_final_distance[:,0]=20
    #ldp_gap_accumulation = np.zeros([len(updated_base_list)+1,len(current_chain)+1])
    seq_gap_accumulation = np.zeros([len(updated_base_list)+1,len(current_chain)+1])
    ldp_gap_accumulation = np.zeros([len(updated_base_list)+1,len(current_chain)+1])
    #init gap count
    #for k in range(1,len(current_chain)+1):
    #    seq_gap_accumulation[0,k]=k
    for k in range(1,len(updated_base_list)+1):
        ldp_gap_accumulation[k,0]=k
    seq_gap_total_limit = int(len(updated_base_list)*0.15)
    ldp_prev_connect = np.zeros([len(updated_base_list)+1,len(current_chain)+1])-1
    ldp_prev_connect[0]=0# which keep calculating reasonable penalty to to avoid open gap.
    # for k in range(1,len(updated_base_list)+1):
    #     ldp_prev_connect[k,:]=k-1
    scratch,pointer,seq_gap_count,ldp_prev_connect,ldp_gap_count=global_align_score(ldp_gap_penalty,seq_gap_penalty, match_matrix,
    np.arange(len(updated_base_list)), np.arange(len(current_chain)),
     pointer, scratch,seq_gap_accumulation, ldp_gap_accumulation,fragment_final_distance, ldp_prev_connect,seq_gap_total_limit)
    if save_path is not None:
        score_path = os.path.join(save_path,"optimal_score.txt")
        np.savetxt(score_path,scratch)
        score_path = os.path.join(save_path,"optimal_direction.txt")
        np.savetxt(score_path,pointer)#verfied corret
        score_path = os.path.join(save_path,"sequence_miss_count.txt")
        np.savetxt(score_path,seq_gap_count)#verified correct
        score_path = os.path.join(save_path,"ldp_prev_connect.txt")
        np.savetxt(score_path,ldp_prev_connect)
        score_path = os.path.join(save_path,"ldp_distance_array.txt")
        np.savetxt(score_path,fragment_final_distance)
        score_path = os.path.join(save_path,"ldp_gap_count.txt")
        np.savetxt(score_path,ldp_gap_count)
    input_seq_line=""
    for i in range(len(updated_base_list)):
        choice=int(np.argmax(updated_base_list[i]))
        input_seq_line += map_dict[choice]
    #trace back all possible choices that with score>0
    #any score<0, we assume directly using the predictions
    max_score_candidate_index = np.argwhere(scratch[len(updated_base_list)]>-1000)
    max_score_list = []
    match_seq_list = []
    match_interval_list =[]
    count_match=0
    for candidate_index in max_score_candidate_index:
        match_matrix = np.zeros(len(updated_base_list))-1
        candidate_index =int(candidate_index)
        max_index = candidate_index
        max_index = [len(updated_base_list),max_index]
        check_pointer = pointer[max_index[0],max_index[1]]
        cur_x=int(max_index[0])
        cur_y=int(max_index[1])
        while check_pointer!=0:
            if check_pointer==999999:
                match_matrix[cur_x-1]=cur_y-1
                current_begin_index = cur_y-1
                cur_x -=1
                cur_y -=1
            elif check_pointer>0:
                cur_y -= int(check_pointer)
            elif check_pointer<0:
                cur_x -= int(abs(check_pointer))
            #match_matrix[cur_x-1]=cur_y-1
            check_pointer = pointer[cur_x,cur_y]
        if save_path is not None:
            score_path = os.path.join(save_path,"match_result%d.txt"%count_match)
            np.savetxt(score_path,match_matrix)
        match_seq_line = ""
        count_replace =0
        for i in range(len(match_matrix)):
            if match_matrix[i]==-1:
                match_seq_line+="-"
            else:
                point_index = int(match_matrix[i])
                #important rule, if the new assignment prob>average, accept, otherwise, decline
                match_label = int(current_chain[point_index])
                current_ldp_score = updated_base_list[i]
                current_ldp_prob = current_ldp_score[match_label]

                match_seq_line+=map_dict[match_label]
                if current_ldp_prob<=0.25:
                    count_replace+=1

        current_score = float(scratch[len(updated_base_list),candidate_index])
        max_score_list.append(current_score)
        match_seq_list.append(match_seq_line)
        match_interval_list.append([current_begin_index,candidate_index])
        count_match+=1
    all_seq_line = ""
    for i in range(len(current_chain)):
        match_label = int(current_chain[i])
        all_seq_line += map_dict[match_label]
    if save_path is not None:
        score_path = os.path.join(save_path,"match_seq.txt")
        with open(score_path,'w') as file:
            file.write(input_seq_line+"\n")
            for j in range(len(match_seq_list)):
                current_interval = match_interval_list[j]
                file.write(match_seq_list[j]+"\t"+"%.2f"%max_score_list[j]+
                "\t%d,%d\n"%(current_interval[0],current_interval[1]))
            file.write(all_seq_line+"\n")

    return max_score_list,match_seq_list,match_interval_list

def greedy_assign_PS(All_Base_Path_List_sugar,All_Path_List_sugar,Path_P_align_list,Path_P_reverse_align_list,
        Pho_Prob_Refer_Dict,Pho_Prob_Refer_Reverse_Dict,chain_prob,
     ldp_size,sugar_point,map_info_list,chain_dict,greedy_save_path,top_select,checking_stride):
    cur_frag_path = os.path.join(greedy_save_path,"All_fragment_%d_%d_%d.pkl"%(ldp_size,checking_stride,top_select))
    cur_frag_path2 = os.path.join(greedy_save_path,"Location_fragment_%d_%d_%d.pkl"%(ldp_size,checking_stride,top_select))
    if os.path.exists(cur_frag_path) and os.path.exists(cur_frag_path2):
        with open(cur_frag_path, 'rb') as handle:
            overall_dict= pickle.load(handle)
        with open(cur_frag_path2, 'rb') as handle:
            frag_location_dict = pickle.load(handle)
        return overall_dict,frag_location_dict
    base_region_prob = chain_prob[2:6]
    merged_cd_dens = sugar_point.merged_cd_dens

    ldp_proper_range = [4,10]
    overall_score=0
    final_check_dir_list=[]
    count_all_combination=0
    count_use_combination=0
    ldp_gap_penalty = 10 # only applied if not reasonable skip
    seq_gap_penalty = 25 # applied once we need to skip an amino acid
    overall_dict=defaultdict(list)#use path id staring index to design overlap
    frag_location_dict = {}

    for path_id,cur_path_list in enumerate(All_Base_Path_List_sugar):
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
        #calculate the reverse score
        #updated_base_score_reverse_list = [updated_base_list[j]+pho_reverse_base_prob_list[j] for j in range(len(current_base_list))]
        #updated_base_score_reverse_list = np.array(updated_base_score_reverse_list)

        cur_path_save_path = os.path.join(greedy_save_path,"path_%d"%path_id)
        mkdir(cur_path_save_path)

        start_ldp_index =0
        count_path_frag = 0
        while start_ldp_index<len(updated_base_list)-5:
            end_ldp_index=min(start_ldp_index+ldp_size,len(updated_base_list))
            study_ldp_base_list = updated_base_score_list[start_ldp_index:end_ldp_index]
            #study_ldp_base_list_reverse = updated_base_score_reverse_list[start_ldp_index:end_ldp_index]
            current_location_list = [merged_cd_dens[int(kk)] for kk in current_path[start_ldp_index:end_ldp_index]]
            fragment_ldp_location= Convert_LDPcoord_To_Reallocation(current_location_list, map_info_list)

            fragment_distance_array = Calculate_Distance_array(fragment_ldp_location)
            path_assign_collection_list=[]
            path_assign_score_list =[]
            path_assign_interval_list =[]
            path_chain_list=[]
            path_direction_list=[]
            path_chainlength_list=[]
            for current_chain_id in chain_dict:
                current_chain = chain_dict[current_chain_id]
                current_reverse_chain = current_chain[::-1]#to match from end to begining condition
                dp_save_path = os.path.join(cur_path_save_path,"path_starting_%d_chain_%s"%(start_ldp_index,current_chain_id))
                mkdir(dp_save_path)
                max_score, match_seq_line,match_seq_interval = dynamic_assign_multi(study_ldp_base_list,current_chain,
                        ldp_gap_penalty,seq_gap_penalty,fragment_distance_array, None)
                path_assign_collection_list.extend(match_seq_line)
                path_assign_score_list.extend(max_score)
                path_assign_interval_list.extend(match_seq_interval)
                path_chain_list.extend([current_chain_id]*len(max_score))
                path_direction_list.extend([1]*len(max_score))
                path_chainlength_list.extend([len(current_chain)]*len(max_score))
                dp_save_path = os.path.join(cur_path_save_path,"rpath_starting_%d_chain_%s"%(start_ldp_index,current_chain_id))
                mkdir(dp_save_path)
                max_score, match_seq_line,match_seq_interval = dynamic_assign_multi(study_ldp_base_list,current_reverse_chain,
                        ldp_gap_penalty,seq_gap_penalty,fragment_distance_array, None)
                path_assign_collection_list.extend(match_seq_line)
                path_assign_score_list.extend(max_score)
                path_assign_interval_list.extend(match_seq_interval)
                path_chain_list.extend([current_chain_id]*len(max_score))
                path_direction_list.extend([-1]*len(max_score))
                path_chainlength_list.extend([len(current_chain)]*len(max_score))

            count_all_combination+=len(path_direction_list)

            path_save_path = os.path.join(cur_path_save_path,"Collecttop_%d_%d"%(path_id,start_ldp_index))
            mkdir(path_save_path)
            if len(os.listdir(path_save_path))>0:
                os.system("rm "+path_save_path+"/*")
            #clear this directory to avoid previous results
            final_check_dir_list.append("Collecttop_%d_%d"%(path_id,start_ldp_index))
            path_assign_score_list=np.array(path_assign_score_list)
            path_score_index = np.argsort(path_assign_score_list)
            path_score_index = path_score_index[::-1]#from bigger to smaller
            rank_id = 1
            final_path_assign_collection_list =[]
            final_path_assign_score_list =[]
            final_path_assign_interval_list =[]
            final_path_chain_list=[]
            final_path_direction_list=[]
            final_path_chainlength_list=[]
            for select_index in path_score_index[:top_select]:
                print("for fragment starting %d, select %d as one of top"%(count_path_frag,select_index))
                current_assign_match_seq = path_assign_collection_list[select_index]
                Visualize_assign_DPbase(path_save_path,"fragement_%d"%rank_id,
                    current_path[start_ldp_index:end_ldp_index],current_assign_match_seq,
                    current_base_list[start_ldp_index:end_ldp_index], sugar_point.merged_cd_dens,map_info_list)
                final_path_assign_collection_list.append(path_assign_collection_list[select_index])
                final_path_assign_score_list.append(path_assign_score_list[select_index])
                final_path_assign_interval_list.append(path_assign_interval_list[select_index])
                final_path_chain_list.append(path_chain_list[select_index])
                final_path_direction_list.append(path_direction_list[select_index])
                final_path_chainlength_list.append(path_chainlength_list[select_index])
                rank_id+=1
            count_use_combination+=min(top_select,len(path_score_index))
            for kk in range(len(final_path_assign_collection_list)):
                new_dict={}
                new_dict['match_seq']=final_path_assign_collection_list[kk]
                new_dict['score']=final_path_assign_score_list[kk]
                new_dict['interval']=final_path_assign_interval_list[kk]
                new_dict['chain']=final_path_chain_list[kk]
                new_dict['direction']=final_path_direction_list[kk]
                # new_dict['ldp_location'] = fragment_ldp_location[kk]
                new_dict['chain_length'] = final_path_chainlength_list[kk]
                overall_dict["%d_%d"%(path_id,start_ldp_index)].append(new_dict)
            frag_location_dict["%d_%d"%(path_id,start_ldp_index)]=fragment_ldp_location
            start_ldp_index+=checking_stride#checking_stride#int(frag_size*0.2)#checking_stride
            count_path_frag+=1
            if len(path_assign_score_list)>0:
                overall_score += np.max(path_assign_score_list)
    print("in total we selected %d/%d possible combinations"%(count_use_combination, count_all_combination))
    with open(cur_frag_path, 'wb') as handle:
        pickle.dump(overall_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    with open(cur_frag_path2, 'wb') as handle:
        pickle.dump(frag_location_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return overall_dict,frag_location_dict

def distributed_dp_calcu(chain_dict,cur_path_save_path,start_ldp_index,
    study_ldp_base_list,ldp_gap_penalty,seq_gap_penalty,fragment_distance_array,
    path_id,top_select):
    path_assign_collection_list=[]
    path_assign_score_list =[]
    path_assign_interval_list =[]
    path_chain_list=[]
    path_direction_list=[]
    path_chainlength_list=[]
    for current_chain_id in chain_dict:
        current_chain = chain_dict[current_chain_id]
        current_reverse_chain = current_chain[::-1]#to match from end to begining condition
        dp_save_path = os.path.join(cur_path_save_path,"path_starting_%d_chain_%s"%(start_ldp_index,current_chain_id))
        mkdir(dp_save_path)
        max_score, match_seq_line,match_seq_interval = dynamic_assign_multi(study_ldp_base_list,current_chain,
                ldp_gap_penalty,seq_gap_penalty,fragment_distance_array, None)
        path_assign_collection_list.extend(match_seq_line)
        path_assign_score_list.extend(max_score)
        path_assign_interval_list.extend(match_seq_interval)
        path_chain_list.extend([current_chain_id]*len(max_score))
        path_direction_list.extend([1]*len(max_score))
        path_chainlength_list.extend([len(current_chain)]*len(max_score))
        #delte all the dirs
        if os.path.exists(dp_save_path):
            os.system("rm -r "+dp_save_path)
        dp_save_path = os.path.join(cur_path_save_path,"rpath_starting_%d_chain_%s"%(start_ldp_index,current_chain_id))
        mkdir(dp_save_path)
        max_score, match_seq_line,match_seq_interval = dynamic_assign_multi(study_ldp_base_list,current_reverse_chain,
                ldp_gap_penalty,seq_gap_penalty,fragment_distance_array, None)
        path_assign_collection_list.extend(match_seq_line)
        path_assign_score_list.extend(max_score)
        path_assign_interval_list.extend(match_seq_interval)
        path_chain_list.extend([current_chain_id]*len(max_score))
        path_direction_list.extend([-1]*len(max_score))
        path_chainlength_list.extend([len(current_chain)]*len(max_score))
        #delte all the dirs
        if os.path.exists(dp_save_path):
            os.system("rm -r "+dp_save_path)

    path_save_path = os.path.join(cur_path_save_path,"Collecttop_%d_%d"%(path_id,start_ldp_index))
    mkdir(path_save_path)
    if len(os.listdir(path_save_path))>0:
        os.system("rm -r "+path_save_path)
    #clear this directory to avoid previous results

    path_assign_score_list=np.array(path_assign_score_list)
    path_score_index = np.argsort(path_assign_score_list)
    path_score_index = path_score_index[::-1]#from bigger to smaller
    rank_id = 1
    final_path_assign_collection_list =[]
    final_path_assign_score_list =[]
    final_path_assign_interval_list =[]
    final_path_chain_list=[]
    final_path_direction_list=[]
    final_path_chainlength_list=[]
    for select_index in path_score_index[:top_select]:
        current_assign_match_seq = path_assign_collection_list[select_index]
        final_path_assign_collection_list.append(path_assign_collection_list[select_index])
        final_path_assign_score_list.append(path_assign_score_list[select_index])
        final_path_assign_interval_list.append(path_assign_interval_list[select_index])
        final_path_chain_list.append(path_chain_list[select_index])
        final_path_direction_list.append(path_direction_list[select_index])
        final_path_chainlength_list.append(path_chainlength_list[select_index])
        rank_id+=1
    if os.path.exists(path_save_path):
        os.system("rm -r "+path_save_path)
    return start_ldp_index,final_path_assign_collection_list,final_path_assign_score_list,\
    final_path_assign_interval_list,final_path_chain_list,final_path_direction_list,final_path_chainlength_list


def greedy_assign_PS_effective(All_Base_Path_List_sugar,All_Path_List_sugar,Path_P_align_list,Path_P_reverse_align_list,
        Pho_Prob_Refer_Dict,Pho_Prob_Refer_Reverse_Dict,chain_prob,
     ldp_size,sugar_point,map_info_list,chain_dict,greedy_save_path,top_select,checking_stride,num_cpus=128):

    cur_frag_path = os.path.join(greedy_save_path,"All_fragment_%d_%d_%d.pkl"%(ldp_size,checking_stride,top_select))
    cur_frag_path2 = os.path.join(greedy_save_path,"Location_fragment_%d_%d_%d.pkl"%(ldp_size,checking_stride,top_select))
    if os.path.exists(cur_frag_path) and os.path.exists(cur_frag_path2):
        with open(cur_frag_path, 'rb') as handle:
            overall_dict= pickle.load(handle)
        with open(cur_frag_path2, 'rb') as handle:
            frag_location_dict = pickle.load(handle)
        return overall_dict,frag_location_dict
    base_region_prob = chain_prob[2:6]
    merged_cd_dens = sugar_point.merged_cd_dens

    ldp_proper_range = [4,10]
    overall_score=0
    #final_check_dir_list=[]
    count_all_combination=0
    count_use_combination=0
    ldp_gap_penalty = 10 # only applied if not reasonable skip
    seq_gap_penalty = 25 # applied once we need to skip an amino acid
    overall_dict=defaultdict(list)#use path id staring index to design overlap
    frag_location_dict = {}


    for path_id,cur_path_list in enumerate(All_Base_Path_List_sugar):
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
        #calculate the reverse score
        #updated_base_score_reverse_list = [updated_base_list[j]+pho_reverse_base_prob_list[j] for j in range(len(current_base_list))]
        #updated_base_score_reverse_list = np.array(updated_base_score_reverse_list)

        cur_path_save_path = os.path.join(greedy_save_path,"path_%d"%path_id)
        mkdir(cur_path_save_path)
        from multiprocessing import Pool
        p= Pool(num_cpus)

        start_ldp_index =0
        count_path_frag = 0
        Res_List=[]
        while start_ldp_index<len(updated_base_list)-5:
            end_ldp_index=min(start_ldp_index+ldp_size,len(updated_base_list))
            study_ldp_base_list = updated_base_score_list[start_ldp_index:end_ldp_index]
            #study_ldp_base_list_reverse = updated_base_score_reverse_list[start_ldp_index:end_ldp_index]
            current_location_list = [merged_cd_dens[int(kk)] for kk in current_path[start_ldp_index:end_ldp_index]]
            fragment_ldp_location= Convert_LDPcoord_To_Reallocation(current_location_list, map_info_list)

            fragment_distance_array = Calculate_Distance_array(fragment_ldp_location)
            res=p.apply_async(distributed_dp_calcu,args=(chain_dict,cur_path_save_path,start_ldp_index,
                study_ldp_base_list,ldp_gap_penalty,seq_gap_penalty,fragment_distance_array,
                path_id,top_select,))
            Res_List.append(res)
            frag_location_dict["%d_%d"%(path_id,start_ldp_index)]=fragment_ldp_location
            start_ldp_index+=checking_stride#checking_stride#int(frag_size*0.2)#checking_stride
            count_path_frag+=1
        p.close()
        p.join()
#get the calculated assignment
        for k in range(len(Res_List)):
            all_info = Res_List[k].get()
            start_ldp_index,final_path_assign_collection_list,final_path_assign_score_list,\
            final_path_assign_interval_list,final_path_chain_list,final_path_direction_list,final_path_chainlength_list=all_info

            for kk in range(len(final_path_assign_collection_list)):
                new_dict={}
                new_dict['match_seq']=final_path_assign_collection_list[kk]
                new_dict['score']=final_path_assign_score_list[kk]
                new_dict['interval']=final_path_assign_interval_list[kk]
                new_dict['chain']=final_path_chain_list[kk]
                new_dict['direction']=final_path_direction_list[kk]
                # new_dict['ldp_location'] = fragment_ldp_location[kk]
                new_dict['chain_length'] = final_path_chainlength_list[kk]
                overall_dict["%d_%d"%(path_id,start_ldp_index)].append(new_dict)
        print("we have %d overall dict"%(len(overall_dict)))


    #print("in total we selected %d/%d possible combinations"%(count_use_combination, count_all_combination))
    with open(cur_frag_path, 'wb') as handle:
        pickle.dump(overall_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    with open(cur_frag_path2, 'wb') as handle:
        pickle.dump(frag_location_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return overall_dict,frag_location_dict
