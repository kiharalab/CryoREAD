from collections import defaultdict
import numpy as np
import os
from graph.LDP_ops import Convert_LDPcoord_To_Reallocation
from ops.os_operation import mkdir
from scipy.spatial.distance import cdist
from ops.math_calcuation import calculate_distance
import random



def build_clusters(all_starting_index_keys,ldp_size):
    all_starting_index_keys.sort()
    cluster_dict=defaultdict(list)
    cluster_id=0
    prev_index = all_starting_index_keys[0]
    cluster_dict[cluster_id].append(prev_index)
    for k in range(1,len(all_starting_index_keys)):
        current_index = all_starting_index_keys[k]
        if prev_index+ldp_size>current_index:
            cluster_dict[cluster_id].append(current_index)
        else:
            cluster_id+=1
            cluster_dict[cluster_id].append(current_index)
        prev_index = current_index
    return cluster_dict

def build_clusters_advanced(all_starting_index_keys,all_seq_indexes):
    all_starting_index_keys.sort()
    cluster_dict=defaultdict(list)
    cluster_id=0
    prev_index = all_starting_index_keys[0]
    cluster_dict[cluster_id].append(prev_index)
    for k in range(1,len(all_starting_index_keys)):
        current_index = all_starting_index_keys[k]
        prev_seq_length = all_seq_indexes[k-1]
        if prev_index+prev_seq_length>current_index:
            cluster_dict[cluster_id].append(current_index)
        else:
            cluster_id+=1
            cluster_dict[cluster_id].append(current_index)
        prev_index = current_index
    return cluster_dict
def build_noclusters_extra(all_starting_index_keys,all_seq_indexes):
    cluster_id=0
    cluster_dict=defaultdict(list)
    for k in range(len(all_starting_index_keys)):
        current_index = all_starting_index_keys[k]
        cluster_dict[cluster_id].append(current_index)
        cluster_id+=1
    return cluster_dict
def verify_chain_check(current_path_dict,current_cluster_id_list):
    current_seq_info = current_path_dict[current_cluster_id_list[0]]
    chain_id = current_seq_info['chain']
    cur_direction = current_seq_info['direction']
    cur_chain_length = current_seq_info['chain_length']
    for k in range(1,len(current_cluster_id_list)):
        current_seq_info = current_path_dict[current_cluster_id_list[k]]
        if current_seq_info['chain']!=chain_id:
            return False
        if current_seq_info['direction']!=cur_direction:
            return False
        if current_seq_info['chain_length']!=cur_chain_length:
            return False
    return True

def merge_cluster_cif_PS(current_path_dict,current_cluster_id_list,ldp_sugar_location,tmp_verify_path,
        Path_P_align_list,Path_P_reverse_align_list,pho_point,map_info_list):
    pho_merged_cd = pho_point.merged_cd_dens[:,:3]
    if len(current_cluster_id_list)==1:
        current_seq_info = current_path_dict[current_cluster_id_list[0]]
        starting_index = int(current_cluster_id_list[0])
        match_seq = current_seq_info['match_seq']
        chain_id = current_seq_info['chain']
        cur_score = current_seq_info['score']
        align_length= len([x for x in match_seq if x!="-"])
        avg_score = cur_score/align_length
        cur_interval = current_seq_info['interval']
        cur_direction = current_seq_info['direction']
        cur_chain_length = current_seq_info['chain_length']
        cur_ldp_size = len(match_seq)
        current_ldp_location = ldp_sugar_location[starting_index:starting_index+cur_ldp_size]
        current_pho_ldp_location = Path_P_align_list[starting_index:starting_index+cur_ldp_size]
        if cur_direction==-1:
            #adjust interval and write from end to begin
            cur_interval = [cur_chain_length-cur_interval[1]-1,cur_chain_length-cur_interval[0]]#reverse chain to correct direction
            match_seq = match_seq[::-1]
            current_ldp_location = current_ldp_location[::-1]
            current_pho_ldp_location = Path_P_reverse_align_list[starting_index:starting_index+cur_ldp_size]
            current_pho_ldp_location = current_pho_ldp_location[::-1]
        #now current_pho_ldp_location are only id list
        current_pho_ldp_align_info = current_pho_ldp_location
        current_pho_ldp_location  = [pho_merged_cd[int(kk)] for kk in current_pho_ldp_location ]
        current_pho_ldp_location = Convert_LDPcoord_To_Reallocation(current_pho_ldp_location , map_info_list)
        #start to ensemble all the information
        fragment_info_list = []
        current_seq_index = cur_interval[0]
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
    print("starting merge fragement process")
    current_cluster_id_list.sort()
    begin_index = current_cluster_id_list[0]
    #cur_ldp_size = len(current_path_dict[current_cluster_id_list[-1]]['match_seq'])
    #end_index = current_cluster_id_list[-1]+cur_ldp_size
    #check all the matches to get the final extended region
    end_index=-1
    for k in range(len(current_cluster_id_list)):
        cur_ldp_size = len(current_path_dict[current_cluster_id_list[k]]['match_seq'])
        cur_end_index = current_cluster_id_list[k]+cur_ldp_size
        if end_index<=cur_end_index:
            end_index=cur_end_index

    overall_range=[begin_index,end_index]
    print("current range:",overall_range)
    assert verify_chain_check(current_path_dict,current_cluster_id_list)#assure same chain, same direction, same chain_length
    current_seq_info = current_path_dict[current_cluster_id_list[0]]
    overall_direction = current_seq_info['direction']
    overall_chain = current_seq_info['chain']
    overall_chain_length = current_seq_info['chain_length']
    #1st fill all the match seq into the begin_index to ending index
    #the rule is use the higher score fragment to avoid double assignment
    overall_match_seq = " "*(end_index-begin_index)
    overall_score_list = [0]*(end_index-begin_index)
    overall_residue_id_list = [-1]*(end_index-begin_index)
    #build a score dict first
    score_dict ={}
    score_list=[]
    for path_id in current_cluster_id_list:
        current_score = current_path_dict[path_id]['score']
        #add a very small random variable to avoid link to same frag
        rand_score = random.random()*1e-4
        score_dict["%.10f"%(current_score+rand_score)]=path_id #starting index
        score_list.append(current_score+rand_score)
    #write low score first, then high score automatically overwrite
    score_sort = np.argsort(score_list)

    for k in range(len(score_sort)):
        sort_index = int(score_sort[k])
        cur_score = score_list[sort_index]
        current_path_id = score_dict["%.10f"%cur_score]
        shift_pos = int(current_path_id-begin_index)
        cur_match_seq = current_path_dict[current_path_id]['match_seq']
        cur_ldp_size = len(cur_match_seq)#in case it's short for some near tail regions
        overall_match_seq=overall_match_seq[:shift_pos]+cur_match_seq+overall_match_seq[(shift_pos+cur_ldp_size):]
        cur_avg_score = cur_score/len([x for x in cur_match_seq if x!="-"])
        for x in range(shift_pos,shift_pos+cur_ldp_size):
            overall_score_list[x]=cur_avg_score
        cur_interval = current_path_dict[current_path_id]['interval']
        current_seq_index = cur_interval[0]
        for kk in range(len(cur_match_seq)):
            if cur_match_seq[kk]=="-":
                continue
            else:
                overall_residue_id_list[shift_pos+kk]=current_seq_index
                current_seq_index+=1

    current_ldp_location = ldp_sugar_location[begin_index:end_index]
    current_pho_ldp_location =  Path_P_align_list[begin_index:end_index]

    if overall_direction==-1:
        overall_match_seq=overall_match_seq[::-1]
        current_ldp_location = current_ldp_location[::-1]
        overall_residue_id_list =overall_residue_id_list[::-1]
        overall_score_list = overall_score_list[::-1]
        #then change the residue id
        overall_residue_id_list=[overall_chain_length-x-1  if x!=-1 else -1 for x in overall_residue_id_list]
        current_pho_ldp_location = Path_P_reverse_align_list[begin_index:end_index]
        current_pho_ldp_location = current_pho_ldp_location[::-1]
    current_pho_ldp_align_info = current_pho_ldp_location
    current_pho_ldp_location  = [pho_merged_cd[int(kk)] for kk in current_pho_ldp_location ]
    current_pho_ldp_location = Convert_LDPcoord_To_Reallocation(current_pho_ldp_location , map_info_list)
    fragment_info_list = []
    overall_score=0
    for k in range(len(overall_match_seq)):
        if overall_match_seq[k]=="-":
            continue
        else:
            cur_location = current_ldp_location[k]
            cur_pho_location = current_pho_ldp_location[k]
            cur_score = overall_score_list[k]
            cur_residue = overall_residue_id_list[k]+1#the original starting from 0, not good
            #make sure assignment in ['A','U','C','G']
            assert overall_match_seq[k] in ['A','U','C','G','T']
            overall_score +=cur_score
            assert cur_residue!=-1
            if  current_pho_ldp_align_info[k]>=0:#avoid non-assign regions always give a wrong point
                fragment_info_list.append([overall_chain,cur_residue,"P",cur_pho_location,overall_match_seq[k],cur_score])
            fragment_info_list.append([overall_chain,cur_residue,"C4'",cur_location,overall_match_seq[k],cur_score])

    #write record for overall output, overall residue id
    #Purpose: Verify our fragment combination is correct
    line_list =[""]*len(overall_match_seq)
    #first write the overall residue type and overall
    for k in range(len(overall_match_seq)):
        line = line_list[k]
        line += overall_match_seq[k]+"\t"
        line += str(overall_residue_id_list[k]+1)+"\t"
        line += "%.2f\t"%overall_score_list[k]
        line_list[k]=line
    #iteratively write different seq
    for path_id in current_cluster_id_list:
        starting_index = int(path_id)
        shift_pos = int(starting_index-begin_index)
        current_seq_info = current_path_dict[path_id]
        match_seq = current_seq_info['match_seq']
        cur_ldp_size = len(match_seq)
        cur_useful_size = len([x for x in match_seq if x!="-"])
        cur_interval = current_seq_info['interval']
        cur_direction = current_seq_info['direction']
        cur_chain_length = current_seq_info['chain_length']
        if cur_direction==-1:
            #adjust interval and write from end to begin
            cur_interval = [cur_chain_length-cur_interval[0]-cur_useful_size,cur_chain_length-cur_interval[0]]#reverse chain to correct direction
            match_seq = match_seq[::-1]
            shift_pos = end_index-begin_index-shift_pos-cur_ldp_size#prev it's the distance to the ending part, now we swapped
        current_seq_index = cur_interval[0]
        cur_score = current_seq_info['score']
        align_length= len([x for x in match_seq if x!="-"])
        avg_score = cur_score/align_length
        for k in range(len(overall_match_seq)):
            line = line_list[k]
            if k>=shift_pos and k<shift_pos+cur_ldp_size:
                if match_seq[k-shift_pos]!="-":
                    line+=match_seq[k-shift_pos]+"\t"
                    line+=str(current_seq_index+1)+"\t"
                    line+="%.2f\t"%(avg_score)
                    current_seq_index+=1
                else:
                    line+="-\t-\t-\t"
            else:
                line+="-\t-\t-\t"
            line_list[k]=line
    with open(tmp_verify_path,'w') as file:
        file.write("Chain Direction:%d\n"%cur_direction)
        for line in line_list:
            file.write(line+"\n")
    return fragment_info_list,overall_score

def clean_pho_assign_info_list(frag_info_list,pho_coordinate,map_info_list,cutoff=7,round=1):
    #clean frag seq_index based on the begining one
    # updated_frag_info_list=[]

    # prev_pho=True
    # for i in range(len(frag_info_list)):
    #     if i==0:
    #         chain_id,current_seq_index,atom_name, cur_atom_position,nuc_type,avg_score =frag_info_list[i]
    #         track_seq_index= current_seq_index
    #         updated_frag_info_list.append(frag_info_list[i])
    #     else:
    #         chain_id,current_seq_index,atom_name, cur_atom_position,nuc_type,avg_score =frag_info_list[i]
    #         if atom_name=="C4'":
    #             if not prev_pho:
    #                 track_seq_index+=1
    #             updated_frag_info_list.append([chain_id,track_seq_index,atom_name, cur_atom_position,nuc_type,avg_score ])
    #             prev_pho=False
    #         else:
    #             track_seq_index+=1
    #             updated_frag_info_list.append([chain_id,track_seq_index,atom_name, cur_atom_position,nuc_type,avg_score ])
    #             prev_pho=True
    # frag_info_list = updated_frag_info_list
    #for duplicated same position, multiple duplications

    Seq_Index_Set={}
    Seq_Index_index_list=defaultdict(list)
    Seq_Index_score_list=defaultdict(list)

    for i in range(len(frag_info_list)):
        chain_id,current_seq_index,atom_name, cur_atom_position,nuc_type,avg_score =frag_info_list[i]
        if current_seq_index not in Seq_Index_Set:
            Seq_Index_Set[current_seq_index]=1
            Seq_Index_index_list[current_seq_index].append(i)
            Seq_Index_score_list[current_seq_index].append(avg_score)
        else:
            Seq_Index_Set[current_seq_index]+=1
            Seq_Index_index_list[current_seq_index].append(i)
            Seq_Index_score_list[current_seq_index].append(avg_score)
    #handling duplicated index error when residue between x-999-y,use different flags to indicate
    error_index_start = 999
    seq_index_list = list(Seq_Index_Set.keys())
    min_seq_index = min(seq_index_list)
    max_seq_index = max(seq_index_list)
    if min_seq_index<=error_index_start and error_index_start<=max_seq_index:
        error_index_start=max_seq_index+1
    if abs(min_seq_index-error_index_start)<=len(Seq_Index_Set):
         error_index_start=max_seq_index+1
    updated_frag_info_list=[]
    for tmp_index in Seq_Index_Set:
        current_index_list = Seq_Index_index_list[tmp_index]
        if len(current_index_list)==1:
            updated_frag_info_list.append(frag_info_list[int(current_index_list[0])])
        else:
            cur_score_list = Seq_Index_score_list[tmp_index]
            max_score = max(cur_score_list)
            unique_score_list = np.unique(cur_score_list)
            sorted_score_index = np.argsort(unique_score_list)
            score_map_index={}
            for k,tmp_score_index in enumerate(sorted_score_index):
                tmp_score = unique_score_list[int(tmp_score_index)]
                score_map_index["%.4f"%tmp_score]= k+1
            for current_check_index in current_index_list:
                current_check_index = int(current_check_index)
                chain_id,current_seq_index,atom_name, cur_atom_position,nuc_type,avg_score =frag_info_list[current_check_index]
                if avg_score==max_score:
                    updated_frag_info_list.append(frag_info_list[current_check_index])
                else:
                    score_decide_index =score_map_index["%.4f"%avg_score]
                    updated_frag_info_list.append([chain_id,error_index_start+score_decide_index,atom_name, cur_atom_position,nuc_type,avg_score])
            error_index_start+=len(score_map_index)-1

    frag_info_list = updated_frag_info_list


    pho_ldp_location = Convert_LDPcoord_To_Reallocation(pho_coordinate , map_info_list)
    Sugar_Record_Dict={}
    Pho_Record_Dict={}
    Pho_Position_List = []
    Sugar_Position_List=[]
    Sugar_IDdefine_dict = {}

    Pho_Seqid_list = []
    for i in range(len(frag_info_list)):
        chain_id,current_seq_index,atom_name, cur_atom_position,nuc_type,avg_score =frag_info_list[i]
        if atom_name=="C4'":
            Sugar_Record_Dict[current_seq_index]=frag_info_list[i]
            Sugar_IDdefine_dict[current_seq_index]=len(Sugar_Position_List)
            Sugar_Position_List.append(cur_atom_position)

        else:
            Pho_Position_List.append(cur_atom_position)
            Pho_Seqid_list.append(current_seq_index)
            Pho_Record_Dict[current_seq_index]= frag_info_list[i]

    Pho_Position_List = np.array(Pho_Position_List)
    Sugar_Position_List = np.array(Sugar_Position_List)
    if len(Pho_Position_List)>0:
        distance_array = cdist(Pho_Position_List,Pho_Position_List)
        Remove_SeqID_set = set()#once in this set,should update pho location to match
        Sugar_ID_list = list(Sugar_Record_Dict.keys())
        remove_pair_list=set()
        for k in range(len(distance_array)):
            current_check_distance = distance_array[k]
            current_seqid = Pho_Seqid_list[k]
            selected_index = np.argwhere(current_check_distance<=0.1)
            for tmp_index in selected_index:
                if tmp_index==k:
                    continue
                next_seqid = Pho_Seqid_list[int(tmp_index)]
                Remove_SeqID_set.add(next_seqid)
                Remove_SeqID_set.add(current_seqid)
                if round==2:
                    min_seq_id = min(current_seqid,next_seqid)
                    max_seq_id = max(current_seqid,next_seqid)
                    remove_pair_list.add("%d_%d"%(min_seq_id,max_seq_id))
    else:
        #no pho assigned, assign now by our closest pho
        distance_array = cdist(Sugar_Position_List,pho_ldp_location)
        final_frag_info_list =[]
        for i in range(len(frag_info_list)):
            chain_id,current_seq_index,atom_name, cur_atom_position,nuc_type,avg_score =frag_info_list[i]
            if atom_name=="C4'":
                final_frag_info_list.append(frag_info_list[i])

            sugar_index= int(Sugar_IDdefine_dict[int(current_seq_index)])
            if sugar_index<1:
                continue
            cur_sp_distance1= distance_array[sugar_index-1]
            cur_sp_distance2= distance_array[sugar_index]
            sp_distance_cur = cur_sp_distance1+cur_sp_distance2
            close_location_index = int(np.argmin(sp_distance_cur))
            if cur_sp_distance1[close_location_index]>cutoff or cur_sp_distance2[close_location_index]>cutoff:
                continue
            select_pho_location = pho_ldp_location[close_location_index]
            cur_info = [chain_id,current_seq_index,'P', select_pho_location,nuc_type,avg_score]
            print("previous loc ",cur_atom_position,"updated loc: ",select_pho_location)
            final_frag_info_list.append(cur_info)
        return final_frag_info_list,True
    print("duplicate use of %d/%d pho positions in atomic"%(len(Remove_SeqID_set),len(Sugar_ID_list)))
    if round==1:
        #check the pho is too far away from any one of 2 consecutive paths
        for i in range(1,len(Sugar_ID_list)):
            current_seq_index = Sugar_ID_list[i-1]
            next_seq_index = Sugar_ID_list[i]
            _,_,_, cur_sugar_position,_,_ =Sugar_Record_Dict[current_seq_index]
            _,_,_, next_sugar_position,_,_ =Sugar_Record_Dict[next_seq_index]
            if next_seq_index not in Pho_Record_Dict:
                Remove_SeqID_set.add(next_seq_index)
                continue
            _,_,_, next_pho_position,_,_ =Pho_Record_Dict[next_seq_index]
            s1_p_distance = calculate_distance(cur_sugar_position,next_pho_position)
            s2_p_distance = calculate_distance(next_sugar_position,next_pho_position)
            if s1_p_distance>=cutoff or s2_p_distance>=cutoff:
                Remove_SeqID_set.add(next_seq_index)
        #recalculated the distance
        print("revise %d/%d pho positions in atomic"%(len(Remove_SeqID_set),len(Sugar_ID_list)))
        print(Remove_SeqID_set)
        distance_array = cdist(Sugar_Position_List,pho_ldp_location)
        final_frag_info_list =[]
        for i in range(len(frag_info_list)):
            chain_id,current_seq_index,atom_name, cur_atom_position,nuc_type,avg_score =frag_info_list[i]
            if atom_name=="C4'":
                final_frag_info_list.append(frag_info_list[i])
                continue

            if current_seq_index not in Remove_SeqID_set:
                final_frag_info_list.append(frag_info_list[i])
                continue
            if i==0:
                final_frag_info_list.append(frag_info_list[i])
                continue
            sugar_index= int(Sugar_IDdefine_dict[int(current_seq_index)])
            cur_sp_distance1= distance_array[sugar_index-1]
            cur_sp_distance2= distance_array[sugar_index]
            sp_distance_cur = cur_sp_distance1+cur_sp_distance2
            close_location_index = int(np.argmin(sp_distance_cur))
            if cur_sp_distance1[close_location_index]>cutoff or cur_sp_distance2[close_location_index]>cutoff:
                continue
            select_pho_location = pho_ldp_location[close_location_index]
            cur_info = [chain_id,current_seq_index,atom_name, select_pho_location,nuc_type,avg_score]
            print("previous loc ",cur_atom_position,"updated loc: ",select_pho_location)
            final_frag_info_list.append(cur_info)

    else:
        #check again removing 1 pho ldp assignment if we still have 2 duplicated pho assignments
        distance_array = cdist(Sugar_Position_List,Pho_Position_List)
        Remove_SeqID_set=set()#pick one out of 2 duplicated phos
        for pair in remove_pair_list:
            remove_seqlist = pair.split("_")
            prev_seqid = int(remove_seqlist[0])
            follow_seqid = int(remove_seqlist[1])

            prev_pho_index = Pho_Seqid_list.index(prev_seqid)
            follow_pho_index = Pho_Seqid_list.index(follow_seqid)
            if prev_pho_index==0:
                #priority to keep the begining pho
                Remove_SeqID_set.add(follow_seqid)
                continue
            sugar_index= int(Sugar_IDdefine_dict[int(prev_seqid)])
            cur_sp_distance1= distance_array[sugar_index-1]
            cur_sp_distance2= distance_array[sugar_index]
            prev_sp_distance = cur_sp_distance1[prev_pho_index]+cur_sp_distance2[prev_pho_index]
            sugar_index= int(Sugar_IDdefine_dict[int(follow_seqid)])
            cur_sp_distance1= distance_array[sugar_index-1]
            cur_sp_distance2= distance_array[sugar_index]
            follow_sp_distance = cur_sp_distance1[follow_pho_index]+cur_sp_distance2[follow_pho_index]

            if prev_sp_distance<=follow_sp_distance:
                Remove_SeqID_set.add(follow_seqid)
            else:
                Remove_SeqID_set.add(prev_seqid)
        print("round 2 we still put %d points into clear set"%(len(Remove_SeqID_set)))
        final_frag_info_list =[]
        for i in range(len(frag_info_list)):
            chain_id,current_seq_index,atom_name, cur_atom_position,nuc_type,avg_score =frag_info_list[i]
            if atom_name=="C4'":
                final_frag_info_list.append(frag_info_list[i])
                continue

            if current_seq_index not in Remove_SeqID_set:
                final_frag_info_list.append(frag_info_list[i])
                continue
        print("clearning changed %d to %d frag info"%(len(frag_info_list),len(final_frag_info_list)))


    return final_frag_info_list,False

