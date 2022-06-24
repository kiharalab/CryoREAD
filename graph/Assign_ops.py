
from scipy.spatial.distance import cdist
from ops.math_calcuation import calculate_distance,calculate_cosine_value
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
