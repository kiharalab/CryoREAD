
from structure.Main_Graph import Main_Graph
from structure.Edge import Edge
import os
from ops.os_operation import mkdir
import numpy as np
import pickle
from scipy.spatial.distance import cdist
from data_processing.map_utils import permute_ns_coord_to_pdb,permute_map_coord_to_pdb
from ops.os_operation import mkdir
from ops.cif_utils import Extract_CIF_coord
from collections import defaultdict
from graph.visualize_utils import visualize_graph
from ops.math_calcuation import calculate_distance

def Pick_Graph_coord(cluster_node_list,merged_cd_dens,graph_edge, Nnode, map_info_list):
    mapc, mapr, maps, origin, nxstart, nystart, nzstart = map_info_list
    #1st filter out all the coordinates in this cluster
    nstart = [nxstart,nystart,nzstart]
    nstart = permute_ns_coord_to_pdb(nstart,mapc,mapr,maps)
    new_origin = [origin[k]+nstart[k] for k in range(3)]
    All_location = []
    ID_dict = {}
    for node_id in cluster_node_list:
        node_id = int(node_id)
        location = merged_cd_dens[node_id,:3]
        location = permute_map_coord_to_pdb(location,mapc,mapr,maps)
        All_location.append([location[k]+new_origin[k] for k in range(3)])
        ID_dict[node_id]= len(All_location)
    Edge_pairs = []


    for i in range(len(graph_edge)):
        id1= graph_edge[i].id1
        id2 = graph_edge[i].id2
        if id1 in cluster_node_list and id2 in cluster_node_list:
            match_id1 = ID_dict[id1]
            match_id2 = ID_dict[id2]
            Edge_pairs.append([match_id1,match_id2])

    print("filtered edge numbers %d"%len(Edge_pairs))
    return All_location,Edge_pairs


def Prune_Graph_Edge(graph,merged_cd_dens,sp_prob,base_prob,prob_threshold=0.1):
    graph_edge=graph.edge
    divide_check =5
    Remain_ID=[]
    for i in range(len(graph_edge)):
        id1= graph_edge[i].id1
        id2 = graph_edge[i].id2
        location1 = merged_cd_dens[id1,:3]
        location2 = merged_cd_dens[id2,:3]
        direction = location2 - location1
        #in the edge, check 5 points to see if they break the rules to decide if we maintain the edge or not
        maintain_label=True
        for k in range(0,divide_check+1):
            current_location = location1+direction*k/divide_check
            x = int(current_location[0])
            y = int(current_location[1])
            z = int(current_location[2])
            #current_sp_prob = sp_prob[x,y,z]
            current_base_prob = base_prob[x,y,z]
            # if current_sp_prob<=prob_threshold:
            #     maintain_label=False
            #     break
            if current_base_prob>=0.5:
                maintain_label=False
                break
        if maintain_label:
            Remain_ID.append(i)
    print("%d/%d edges are remained after prob map checking"%(len(Remain_ID),len(graph_edge)))
    return Remain_ID

def verify_edge_criteria(location1,location2,base_prob,sp_prob,prob_threshold):
    divide_check =5
    direction = location2 - location1
        #in the edge, check 5 points to see if they break the rules to decide if we maintain the edge or not
    maintain_label=True
    for k in range(0,divide_check+1):
        current_location = location1+direction*k/divide_check
        x = int(current_location[0])
        y = int(current_location[1])
        z = int(current_location[2])
        current_sp_prob = sp_prob[x,y,z]
        current_base_prob = base_prob[x,y,z]
        if current_sp_prob<=prob_threshold:
            maintain_label=False
            break
        if current_base_prob>=0.5:
            maintain_label=False
            return False

    return maintain_label

def extend_graph_edge_connect(pho_graph,merged_cd_dens,sp_prob,base_prob,tmp_save_path,cutoff,prob_threshold=0.1):
    node_pair_path = os.path.join(tmp_save_path,"Edge_extend_info.txt")

    if not os.path.exists(node_pair_path):
        distance_array = cdist(merged_cd_dens,merged_cd_dens)
        #get all the connected information
        connect_info=set()
        graph_edge = pho_graph.edge
        for i in range(len(graph_edge)):
            id1= graph_edge[i].id1
            id2 = graph_edge[i].id2
            min_id = min(id1,id2)
            max_id = max(id1,id2)
            connect_info.add("%d_%d"%(min_id,max_id))

        Add_Info=[]#[id1,id2,distance]
        for k in range(len(distance_array)):
            current_distance = distance_array[k]
            connect_index_list = np.argwhere(current_distance<cutoff)
            for tmp_index in connect_index_list:
                #verify if it's already in connection list
                tmp_index = int(tmp_index)
                min_index = min(k,tmp_index)
                max_index = max(k,tmp_index)
                if k==tmp_index:
                    continue
                search_key = "%d_%d"%(min_index,max_index)
                if search_key not in connect_info:
                    #verify if the connection cross the base region
                    maintain_label = verify_edge_criteria(merged_cd_dens[k],merged_cd_dens[tmp_index],base_prob,sp_prob,prob_threshold)
                    if maintain_label:
                        Add_Info.append([k,tmp_index,distance_array[k,tmp_index]])
                        connect_info.add(search_key)
        Add_Info = np.array(Add_Info)
        np.savetxt(node_pair_path,Add_Info)
    else:
        Add_Info = np.loadtxt(node_pair_path)
    #use this to add edges to graph.edge
    print("Further adding %d edges based on extending and checking base region"%len(Add_Info))
    for i in range(len(Add_Info)):
        tmp_edge = Edge()
        tmp_edge.id1= int(Add_Info[i][0])
        tmp_edge.id2= int(Add_Info[i][1])
        tmp_edge.d= Add_Info[i][2]*3#penalize those extended connections
        pho_graph.edge.append(tmp_edge)


def Filter_subgraph_edge_density(cluster_node_list,graph_edge):
    edge_density = []

    for i in range(len(graph_edge)):
        id1= graph_edge[i].id1
        id2 = graph_edge[i].id2
        if id1 in cluster_node_list and id2 in cluster_node_list:
            edge_density.append(graph_edge[i].dens)
    return edge_density


def Trace_Int_Location_Dict(merged_cd_dens):
    search_dict = defaultdict(list) #[int_location]:[node_id] list
    actual_location_dict = defaultdict(list) #[int_location]:[actual_location] list
    for k in range(len(merged_cd_dens)):
        cur_location = merged_cd_dens[k]
        x,y,z = cur_location[0],cur_location[1],cur_location[2]
        int_x,int_y,int_z = int(x),int(y),int(z)
        search_key = str(int_x)+","+str(int_y)+","+str(int_z)
        search_dict[search_key].append(k)
        actual_location_dict[search_key].append([x,y,z])
    return search_dict,actual_location_dict

def Identify_Trace_Node(merged_cd_dens,extract_coord_list,map_info_list):

    Filter_search_dict,Filter_location_dict = Trace_Int_Location_Dict(extract_coord_list)

    mapc, mapr, maps, origin, nxstart, nystart, nzstart = map_info_list
    #1st filter out all the coordinates in this cluster
    nstart = [nxstart,nystart,nzstart]
    nstart = permute_ns_coord_to_pdb(nstart,mapc,mapr,maps)
    new_origin = [origin[k]+nstart[k] for k in range(3)]
    modified_merged_cd = []
    for node_id in range(len(merged_cd_dens)):
        node_id = int(node_id)
        location = merged_cd_dens[node_id,:3]
        location = permute_map_coord_to_pdb(location,mapc,mapr,maps)
        new_location = [float("%8.3f"%(location[k]+new_origin[k])) for k in range(3)]#very important for dict search to keep same format
        modified_merged_cd.append(new_location)

    GT_search_dict,GT_location_dict = Trace_Int_Location_Dict(modified_merged_cd)
    print("Extract search dict %d, GT search dict %d"%(len(Filter_search_dict),len(GT_search_dict)))
    Map_Dict = {}#[path_id]:[node_id]
    for key in Filter_search_dict:
        current_id_list = Filter_search_dict[key]
        current_location_list = Filter_location_dict[key]
        gt_id_list = GT_search_dict[key]
        gt_location_list = GT_location_dict[key]
        #assert  len(current_id_list)==len(gt_id_list)
        for k in range(len(current_id_list)):
            current_id = current_id_list[k]
            current_location = current_location_list[k]
            for j in range(len(gt_id_list)):
                gt_id = gt_id_list[j]
                gt_location = gt_location_list[j]
                distance = calculate_distance(current_location,gt_location)
                if distance<0.01:
                    Map_Dict[current_id]=gt_id
                    break
    print("search length %d, return search result length %d"%(len(extract_coord_list),len(Map_Dict)))
    assert len(extract_coord_list)==len(Map_Dict)
    key_list = list(Map_Dict.keys())
    key_list.sort()
    Final_list =[]
    for key in key_list:
        current_node_id = Map_Dict[key]
        Final_list.append(current_node_id)
    return Final_list
def Filter_Path_With_Edge(coordinate_list,edge_pairs,cut_off_length):

    select_lists = []
    #begin_flag=False
    tmp_list=[coordinate_list[0]]
    #prev_cid = cid[0]
    for k in range(1,len(coordinate_list)):
        #cur_cid = cid[k]
        prev_coord = coordinate_list[k-1]
        cur_coord = coordinate_list[k]
        cur_distance = calculate_distance(prev_coord,cur_coord)
        if cur_distance<cut_off_length:
            tmp_list.append(coordinate_list[k])
        if cur_distance>=cut_off_length:
            select_lists.append(tmp_list)
            tmp_list = [coordinate_list[k]]
    if len(tmp_list)>1:
        select_lists.append(tmp_list)
    return select_lists
def Prune_Selected_Path(tmp_save_path,listfiles,merged_cd_dens,drna_graph,map_info_list,ext_name,cut_off_length):
    split_save_path = os.path.join(tmp_save_path,"Prune")
    mkdir(split_save_path)
    count_iter=0
    All_Node_Path_List=[]
    for j,item in enumerate(listfiles):
        cur_cif_path = os.path.join(tmp_save_path,item)
        #extract the coordinates
        extract_coord_list = Extract_CIF_coord(cur_cif_path) #order indicates connections
        if len(extract_coord_list)<=2:
            continue
        #identify the corresponding node id
        node_id_list = Identify_Trace_Node(merged_cd_dens,extract_coord_list,map_info_list)
        All_Node_Path_List.append(node_id_list)
        #first visualize the selected node
        coordinate_list, edge_pairs = Pick_Graph_coord(node_id_list, merged_cd_dens,
                                               drna_graph.edge, drna_graph.Nnode, map_info_list)

        visualize_graph(split_save_path,ext_name+"_graph%d"%j,coordinate_list,edge_pairs)
        select_path_list = Filter_Path_With_Edge(coordinate_list,edge_pairs,cut_off_length)
        print("%d node path find %d sub path"%(len(coordinate_list),len(select_path_list)))
        for select_path in select_path_list:
            visualize_graph(split_save_path,ext_name+"_path%d"%count_iter,select_path,None)
            count_iter+=1
    return All_Node_Path_List

def construct_graph(input_mrc,input_density,pho_point,sp_prob,base_prob,save_path,ext_name,map_info_list,params,prob_threshold=0.1,extend=True):
    #construct the graph to build edges for further processing.
    pho_graph = Main_Graph(params)
    tmp_save_path = os.path.join(save_path,ext_name)
    mkdir(tmp_save_path)
    pho_graph.setup_graph(pho_point, tmp_save_path)
    coordinate_list, edge_pairs = Pick_Graph_coord(list(np.arange(pho_graph.Nnode)), pho_point.merged_cd_dens,
                                                       pho_graph.edge, pho_graph.Nnode, map_info_list)
    visualize_graph(tmp_save_path,ext_name,coordinate_list,edge_pairs)
    #1.3 prune graph edges by prediction map
    #rule 1: phosphate_prob<=prob_threshold
    #rule 2: do not cross base regions >0.5
    edge_prune_path = os.path.join(tmp_save_path,"edge_remain_id.pkl")
    if not os.path.exists(edge_prune_path):
        Remain_ID_List= Prune_Graph_Edge(pho_graph,pho_point.merged_cd_dens,sp_prob,base_prob,prob_threshold=prob_threshold)
        with open(edge_prune_path, 'wb') as handle:
            pickle.dump(Remain_ID_List, handle, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        with open(edge_prune_path, 'rb') as handle:
            Remain_ID_List= pickle.load(handle)
    pho_graph.edge = [edge for k,edge in enumerate(pho_graph.edge) if k in Remain_ID_List]
    coordinate_list, edge_pairs = Pick_Graph_coord(list(np.arange(pho_graph.Nnode)), pho_point.merged_cd_dens,
                                                       pho_graph.edge, pho_graph.Nnode, map_info_list)
    visualize_graph(tmp_save_path,ext_name+"_prune",coordinate_list,edge_pairs)
    #add more edges by using params['R'] to connect all nodes that within R while did not pass any base_region>0.5 and phosphate_prob>=threshold
    if extend:
        extend_graph_edge_connect(pho_graph,pho_point.merged_cd_dens,sp_prob, base_prob,tmp_save_path,params['R'],prob_threshold)
        coordinate_list, edge_pairs = Pick_Graph_coord(list(np.arange(pho_graph.Nnode)), pho_point.merged_cd_dens,
                                                        pho_graph.edge, pho_graph.Nnode, map_info_list)
        visualize_graph(tmp_save_path,ext_name+"_extend",coordinate_list,edge_pairs)
    pho_graph.Ne = len(pho_graph.edge)
    edge_d_dens= pho_graph.set_edge_dens(tmp_save_path,input_mrc,pho_point,input_density)
    return pho_graph,coordinate_list,edge_pairs,edge_d_dens
