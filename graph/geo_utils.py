
import imp
from scipy.spatial.distance import cdist
import numpy as np
from collections import defaultdict

from graph.LDP_ops import permute_point_coord_to_global_coord



def Match_Sugar_Base_Location(Path_ID_List,sugar_point,Base_LDP_List,map_info_list):
    #a dict that matches [global_sugar_coord]:[global_assigned_base_coord]
    sugar_coordinate = sugar_point.merged_cd_dens[:,:3]

    Base_Coord_List = []
    for base_ldp in Base_LDP_List:
        if len(base_ldp.merged_cd_dens)>0:
            base_coord = base_ldp.merged_cd_dens[:,:3]
            Base_Coord_List.append(base_coord)
    Base_Coord_List = np.concatenate(Base_Coord_List,axis=0)
    sb_distance = cdist(sugar_coordinate,Base_Coord_List)
    Base_refer_dict={}
    for cur_path_list in Path_ID_List:
        current_length = len(cur_path_list)
        for k in range(current_length):
            node_id1 = int(cur_path_list[k])

            node_id1_sugar_location = sugar_coordinate[node_id1]
            check_location = ""
            for kk in range(3):
                check_location+="%.4f,"%node_id1_sugar_location[kk]
            sb_distance_node1 = sb_distance[node_id1]
            nearby_index = np.argmin(sb_distance_node1)
            cur_base_location = Base_Coord_List[nearby_index]
            global_sugar_coord = permute_point_coord_to_global_coord(node_id1_sugar_location,map_info_list)
            new_key = ""
            for k in range(3):
                new_key+="%.4f,"%(global_sugar_coord[k])
            global_base_coord = permute_point_coord_to_global_coord(cur_base_location,map_info_list)
            Base_refer_dict[new_key]=global_base_coord
    return Base_refer_dict
