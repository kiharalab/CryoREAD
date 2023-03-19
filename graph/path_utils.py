from ops.os_operation import mkdir
import pickle
import os
from graph.Graph_ops import Pick_Graph_coord,Filter_subgraph_edge_density,Prune_Selected_Path
from graph.visualize_utils import visualize_graph
from graph.ortool_path_ops import ortools_build_path

def collect_all_searched_path(subgraphs, pho_point,pho_graph,search_dir,map_info_list,params):


    all_path_file = os.path.join(search_dir,"All_path.pkl")
    if os.path.exists(all_path_file):
        with open(all_path_file, 'rb') as handle:
            All_Path_List= pickle.load(handle)
        return All_Path_List

    All_Path_List = []
    for cluster_id, subgraph in enumerate(subgraphs):
        if len(subgraph)<5:
            continue
        subgraph_path = os.path.join(search_dir,"sub_graph%d"%cluster_id)
        mkdir(subgraph_path)
        pho_coordinate_list, pho_edge_pairs = Pick_Graph_coord(subgraph, pho_point.merged_cd_dens,
                                                       pho_graph.edge, pho_graph.Nnode, map_info_list)
        visualize_graph(subgraph_path,"pho_graph%d"%cluster_id,pho_coordinate_list,pho_edge_pairs)
        pho_edge_d_dens = Filter_subgraph_edge_density(subgraph,pho_graph.edge)
        assert len(pho_edge_d_dens)==len(pho_edge_pairs)
        listfiles = [x for x in os.listdir(subgraph_path) if ".cif" in x and 'path' in x]
        listfiles.sort()
        if len(listfiles)<1:
            ortools_build_path(subgraph_path,pho_coordinate_list,pho_edge_pairs,pho_edge_d_dens,params['R'],relax_choice=True)
            listfiles = [x for x in os.listdir(subgraph_path) if ".cif" in x and 'path' in x]
            listfiles.sort()
        Path_ID_List=Prune_Selected_Path(subgraph_path,listfiles,pho_point.merged_cd_dens,
                                pho_graph,map_info_list,"pho_prune",params['R'])
        for item in Path_ID_List:
            All_Path_List.append(item)
    with open(all_path_file, 'wb') as handle:
        pickle.dump(All_Path_List, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return All_Path_List
