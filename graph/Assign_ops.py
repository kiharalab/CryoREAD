
from scipy.spatial.distance import cdist
from ops.math_calcuation import calculate_distance,calculate_cosine_value
import numpy as np
def Assign_PhoLDP_Sugarpath(Path_ID_List,sugar_point,pho_point,cut_off_length=7):
    sugar_coodinate = sugar_point.merged_cd_dens[:,:3]
    pho_coordinate = pho_point.merged_cd_dens[:,:3]
    distance_array = cdist(sugar_coodinate,pho_coordinate)
    Path_P_align_list=[]
    Path_P_reverse_align_list=[]
    Path_Pho_ID_List=[]

    for cur_path_list in Path_ID_List:
        current_length = len(cur_path_list)
        tmp_point_list = np.zeros([current_length,2])-1#use -1 to indicate the non-assignment
        for k in range(current_length-1):
            node_id1 = int(cur_path_list[k])
            node_id2 = int(cur_path_list[k+1])
            distance_node1 = distance_array[node_id1]
            distance_node2 = distance_array[node_id2]
            combine_distance = distance_node1+distance_node2
            nearby_index = int(np.argmin(combine_distance))
            minimum_distance_now = combine_distance[nearby_index]
            if minimum_distance_now>cut_off_length*2:
                continue



            location_node1 = sugar_coodinate[node_id1]
            location_node2 = sugar_coodinate[node_id2]
            location_pho_now = pho_coordinate[nearby_index]

            s1_s2_edge = calculate_distance(location_node1,location_node2)
            s1_p_edge = distance_node1[nearby_index]
            s2_p_edge = distance_node2[nearby_index]

            cosine_s2_s1_p = calculate_cosine_value(s1_s2_edge,s1_p_edge,s2_p_edge)
            cosine_s1_s2_p = calculate_cosine_value(s1_s2_edge,s2_p_edge,s1_p_edge)
            cosine_s1_p_s2 = calculate_cosine_value(s2_p_edge,s1_p_edge,s1_s2_edge)
            if cosine_s2_s1_p<0 or cosine_s1_p_s2<0 or cosine_s1_s2_p<0:
                continue
            #all angles<=90, we can apply projection to the main path
            tmp_point_list[k,1]= nearby_index#end pointer for this
            tmp_point_list[k+1,0]=nearby_index #starting pointer for the following point
        #assign prev pho node for 1st node , following node of the last node
        begining_node = cur_path_list[0]
        distance_begin_node = distance_array[begining_node]
        distance_next_node = distance_array[int(cur_path_list[1])]
        location_node1 = sugar_coodinate[begining_node]
        location_node2 = sugar_coodinate[int(cur_path_list[1])]
        possible_candidate_index = np.argwhere(distance_begin_node<=cut_off_length)
        for candidate_index in possible_candidate_index:
            s1_s2_edge = calculate_distance(location_node1,location_node2)
            s1_p_edge = distance_begin_node[candidate_index]
            s2_p_edge = distance_next_node[candidate_index]
            cosine_s2_s1_p = calculate_cosine_value(s1_s2_edge,s1_p_edge,s2_p_edge)
            if cosine_s2_s1_p<0:
                tmp_point_list[0,0]=candidate_index
                break
        #assign the pho node for the last sugar node in the path
        ending_node = cur_path_list[-1]
        distance_end_node = distance_array[ending_node]
        distance_prev_node = distance_array[int(cur_path_list[-2])]
        location_node1 = sugar_coodinate[ending_node]
        location_node2 = sugar_coodinate[int(cur_path_list[-2])]
        possible_candidate_index = np.argwhere(distance_end_node<=cut_off_length)
        second_choice = -1

        for candidate_index in possible_candidate_index:
            s1_s2_edge = calculate_distance(location_node1,location_node2)
            s1_p_edge = distance_end_node[candidate_index]
            s2_p_edge = distance_prev_node[candidate_index]
            cosine_s2_s1_p = calculate_cosine_value(s1_s2_edge,s1_p_edge,s2_p_edge)
            if cosine_s2_s1_p<0:
                tmp_point_list[-1,1]=candidate_index
                break



        #assign previous pointer to every allowed sugar points as a reward for them
        #from begining to ending, assign starting point
        cur_align_list1=[]
        cur_id_list = []
        cur_starting_point=tmp_point_list[0,0]
        for k in range(len(tmp_point_list)):
            if tmp_point_list[k,0]==-1:
                node_id1 = int(cur_path_list[k])
                distance_node1 = distance_array[node_id1]
                cur_distance = distance_node1[int(cur_starting_point)]
                if cur_distance<=cut_off_length:
                    cur_align_list1.append(cur_starting_point)
                else:
                    cur_align_list1.append(-1)
            else:
                cur_starting_point=tmp_point_list[k,0]
                cur_align_list1.append(cur_starting_point)
                cur_id_list.append(cur_starting_point)
        cur_id_list.append(tmp_point_list[-1,1])
        Path_Pho_ID_List.append(cur_id_list)
        #2 from ending to begining, assign ending point
        cur_reverse_align_list =[]
        cur_starting_point=tmp_point_list[-1,1]
        for k in range(len(tmp_point_list)-1, -1, -1):
            if tmp_point_list[k,1]==-1:
                node_id1 = int(cur_path_list[k])
                distance_node1 = distance_array[node_id1]
                cur_distance = distance_node1[int(cur_starting_point)]
                if cur_distance<=cut_off_length:
                    cur_reverse_align_list.append(cur_starting_point)
                else:
                    cur_reverse_align_list.append(-1)
            else:
                cur_starting_point=tmp_point_list[k,1]
                cur_reverse_align_list.append(cur_starting_point)
        cur_reverse_align_list=cur_reverse_align_list[::-1]
        assert len(cur_align_list1)==len(cur_reverse_align_list) and len(cur_align_list1)==len(cur_path_list)
        Path_P_align_list.append(cur_align_list1)
        Path_P_reverse_align_list.append(cur_reverse_align_list)
    return Path_P_align_list,Path_P_reverse_align_list,Path_Pho_ID_List
def Assign_Base_Main_Path_sugar(Path_ID_List,pho_coordinate,Base_LDP_List,base_prob_array,cut_off_length=10):

    Base_Coord_List = []
    for base_ldp in Base_LDP_List:
        if len(base_ldp.merged_cd_dens)>0:
            base_coord = base_ldp.merged_cd_dens[:,:3]
            Base_Coord_List.append(base_coord)
    Base_Coord_List = np.concatenate(Base_Coord_List,axis=0)

    distance_array = cdist(pho_coordinate,Base_Coord_List)
    All_Base_Assign_List =[]#list of list: where probability of specific base will be put here.
    count_out =0
    count_total=0
    for cur_path_list in Path_ID_List:
        current_length = len(cur_path_list)
        count_total+=current_length
        tmp_list = []
        for k in range(current_length):
            node_id1 = int(cur_path_list[k])
            distance_node1 = distance_array[node_id1]
            nearby_index = np.argmin(distance_node1)
            minimum_distance_now = distance_node1[nearby_index]
            cur_base_coord = Base_Coord_List[int(nearby_index)]
            current_tmp_prob = base_prob_array[:,int(cur_base_coord[0]),int(cur_base_coord[1]),int(cur_base_coord[2])]
            if minimum_distance_now>cut_off_length:
                count_out+=1
                tmp_list.append(current_tmp_prob*0.5)
                continue
            else:
                tmp_list.append(current_tmp_prob)
        All_Base_Assign_List.append(tmp_list)
    print("we have %d/%d base assignment outside safe range: %.3f"%(count_out,count_total,count_out/count_total))
    return All_Base_Assign_List


def Assign_Base_Main_Path(Path_ID_List,pho_coordinate,Base_LDP_List,base_prob_array,cut_off_length=10,reverse_flag=False,return_dict=False):

    Base_Coord_List = []
    for base_ldp in Base_LDP_List:
        if len(base_ldp.merged_cd_dens)>0:
            base_coord = base_ldp.merged_cd_dens[:,:3]
            Base_Coord_List.append(base_coord)
    Base_Coord_List = np.concatenate(Base_Coord_List,axis=0)

    distance_array = cdist(pho_coordinate,Base_Coord_List)
    All_Base_Assign_List =[]#list of list: where probability of specific base will be put here.
    Assign_Refer_Dict = {}#[node_id]:[prob]
    count_out =0
    count_total=0
    for cur_path_list in Path_ID_List:
        current_length = len(cur_path_list)
        count_total+=current_length
        if reverse_flag is True:
            cur_path_list = cur_path_list[::-1]
        tmp_list = []
        for k in range(current_length-1):
            node_id1 = int(cur_path_list[k])
            node_id2 = int(cur_path_list[k+1])
            distance_node1 = distance_array[node_id1]
            distance_node2 = distance_array[node_id2]
            combine_distance = distance_node1+distance_node2
            nearby_index = np.argmin(combine_distance)
            minimum_distance_now = combine_distance[nearby_index]
            cur_base_coord = Base_Coord_List[int(nearby_index)]
            current_tmp_prob = base_prob_array[:,int(cur_base_coord[0]),int(cur_base_coord[1]),int(cur_base_coord[2])]
            if minimum_distance_now>cut_off_length*2:
                count_out+=1
                tmp_list.append(current_tmp_prob*0.5)
                Assign_Refer_Dict[node_id1]=current_tmp_prob*0.5
                continue
            else:
                tmp_list.append(current_tmp_prob)
                Assign_Refer_Dict[node_id1]=current_tmp_prob
        if reverse_flag:
            tmp_list = tmp_list[::-1]
        All_Base_Assign_List.append(tmp_list)
    print("we have %d/%d base assignment outside safe range: %.3f"%(count_out,count_total,count_out/count_total))
    if return_dict:
        return All_Base_Assign_List,Assign_Refer_Dict
    else:
        return All_Base_Assign_List
