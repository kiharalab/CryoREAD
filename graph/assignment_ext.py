import numpy as np
import os

import pickle
from collections import defaultdict
def identify_last_match(match_seq):
    match_index=0
    for k in range(len(match_seq)):
        if match_seq[k]=="-":
            continue
        match_index=k
    return match_index
def reorganize_missing_list(missing_list,define_ldp_size):
    final_missing_list = []
    for k in range(len(missing_list)):
        start,end = missing_list[k]
        if end-start<=define_ldp_size:
            final_missing_list.append([start,end])
        else:
            gap_range=end-start
            # divide_num = int(np.ceil(gap_range/define_ldp_size))
            # small_gap = int(gap_range/divide_num)
            # print("gap design %d, we need to go through %d frags"%(small_gap,divide_num))
            # print("previous grap %d/%d"%(start,end))
            # for kk in range(divide_num):
            #     if kk!=divide_num-1:
            #         final_missing_list.append([start,start+small_gap-1])
            #         start+=small_gap
            #     else:
            #         final_missing_list.append([start,end])
            max_check_length= define_ldp_size-1
            while start<end:
                cur_end = min(start+max_check_length,end)
                final_missing_list.append([start,cur_end])
                start+=max_check_length+1

    print(final_missing_list)
    print("reformed missing list length %d, previous %d"%(len(final_missing_list),len(missing_list)))
    return final_missing_list


def Extend_Solve_Assignment_SP_support(save_path,All_Path_List_sugar,solve_frag_combine_list,
order_key_index,order_chain_index,overall_dict,gap_limit=3,define_ldp_size=20):
    #1 identify all the fragments that has been identified by current methods
    occupy_dict= defaultdict(list)
    Match_Dict = defaultdict(list)
    extend_support_path = os.path.join(save_path,"SP_fillextra_info.pkl")
    if os.path.exists(extend_support_path):
        with open(extend_support_path, 'rb') as handle:
            external_dict= pickle.load(handle)

        return external_dict

    for order_index in solve_frag_combine_list:
        order_index = int(order_index)
        current_key1 = order_key_index[order_index]
        current_chain_candidate1 = order_chain_index[order_index]
        split_key1 = current_key1.split("_")
        ldp_path1= int(split_key1[0])
        ldp_starting_index1 = int(split_key1[1])
        current_seq_info1 = overall_dict[current_key1][current_chain_candidate1]
        occupy_dict[ldp_path1].append([ldp_starting_index1,ldp_starting_index1+len(current_seq_info1['match_seq'])])
        Match_Dict[ldp_path1].append(current_seq_info1)
    #2 iterative check remaining regions without assignment
    Missing_Range_Dict=defaultdict(list)
    for ldp_path_id in range(len(All_Path_List_sugar)):
        cur_sugar_ldp_path_list = All_Path_List_sugar[ldp_path_id]
        current_path_match_info = Match_Dict[ldp_path_id]
        identified_ldp_frag_list = occupy_dict[ldp_path_id]
        if len(identified_ldp_frag_list)==0:
            #case that no frag identified for this path
            if len(cur_sugar_ldp_path_list)<=gap_limit:
                continue
            Missing_Range_Dict[ldp_path_id].append([0,len(cur_sugar_ldp_path_list)])
            continue
        identified_ldp_frag_list = sorted(identified_ldp_frag_list,key=(lambda x:x[0]))
        for k in range(len(identified_ldp_frag_list)-1):
            if k==0:
                starting_index,_=identified_ldp_frag_list[0]
                if starting_index>gap_limit:
                    Missing_Range_Dict[ldp_path_id].append([0,starting_index])
            else:
                cur_starting_index,cur_ending_index = identified_ldp_frag_list[k]
                next_starting_index,next_ending_index = identified_ldp_frag_list[k+1]
                if next_starting_index-cur_ending_index<=gap_limit:
                    continue
                current_match_info = current_path_match_info[k]#check which index is the last matched residue, then we use that as a starting referece
                cur_last_real_match_idx = identify_last_match(current_match_info['match_seq'])+cur_starting_index
                Missing_Range_Dict[ldp_path_id].append([cur_last_real_match_idx,next_starting_index])
        total_ldp_list_length = len(cur_sugar_ldp_path_list)
        final_ending_frag_index = identified_ldp_frag_list[-1][1]
        if total_ldp_list_length-final_ending_frag_index>gap_limit:
            current_match_info = current_path_match_info[-1]#check which index is the last matched residue, then we use that as a starting referece
            cur_last_real_match_idx = identify_last_match(current_match_info['match_seq'])+identified_ldp_frag_list[-1][0]
            Missing_Range_Dict[ldp_path_id].append([cur_last_real_match_idx,total_ldp_list_length])
    Path_Assign_Dict=defaultdict(dict)
    #3 fill the missing alignment region simply use dp results.
    for ldp_path_id in range(len(All_Path_List_sugar)):
        cur_sugar_ldp_path_list = All_Path_List_sugar[ldp_path_id]
        missing_list = Missing_Range_Dict[ldp_path_id]
        missing_list = reorganize_missing_list(missing_list,define_ldp_size)#specifically for strict since some missing region is even bigger compared to
        for k in range(len(missing_list)):
            fill_starting_index, fill_ending_index = missing_list[k]
            current_key = "%d_%d"%(ldp_path_id,int(fill_starting_index))

            if len(cur_sugar_ldp_path_list)-fill_starting_index<=gap_limit:
                continue
            current_seq_info = overall_dict[current_key]
            iter_check=0
            while len(current_seq_info)==0 and iter_check<=5:
                fill_starting_index-=1
                current_key = "%d_%d"%(ldp_path_id,int(fill_starting_index))
                current_seq_info = overall_dict[current_key]
                iter_check+=1
            if iter_check>5:
                #search another direction
                iter_check=0
                while len(current_seq_info)==0 and iter_check<=5:
                    fill_starting_index+=1
                    current_key = "%d_%d"%(ldp_path_id,int(fill_starting_index))
                    current_seq_info = overall_dict[current_key]
                    iter_check+=1
            if iter_check>5:
                continue
            print("extend for %d-%d"%(fill_starting_index,fill_ending_index))
            current_seq_info = overall_dict[current_key][0]#use top 1 as match candidate
            print("before info:",current_seq_info)
            remove_match_seq = current_seq_info['match_seq'][(fill_ending_index-fill_starting_index):]
            current_seq_info['match_seq']=current_seq_info['match_seq'][:fill_ending_index-fill_starting_index]

            count_remove_match_num = len([1 for tmp_match_id in remove_match_seq if tmp_match_id!="-"])

            current_seq_info['interval'][1]=current_seq_info['interval'][1]-count_remove_match_num
            print("previous info:",current_seq_info)
            Path_Assign_Dict[ldp_path_id][fill_starting_index]= current_seq_info
    with open(extend_support_path, 'wb') as handle:
        pickle.dump(Path_Assign_Dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return Path_Assign_Dict





