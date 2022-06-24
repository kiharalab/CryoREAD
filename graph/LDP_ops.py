
import os
import numpy as np
from graph.Points import Points
from data_processing.map_utils import permute_ns_coord_to_pdb,permute_map_coord_to_pdb
from graph.io_utils import save_LDP_map
from graph.visualize_utils import Show_Graph_Connect,Show_Bfactor_cif
def Extract_LDP_coord(merged_cd_dens, mapc, mapr, maps, origin, nxstart, nystart, nzstart):
    #1st filter out all the coordinates in this cluster
    nstart = [nxstart,nystart,nzstart]
    nstart = permute_ns_coord_to_pdb(nstart,mapc,mapr,maps)
    new_origin = [origin[k]+nstart[k] for k in range(3)]
    All_location = []
    for node_id in range(len(merged_cd_dens)):
        node_id = int(node_id)
        location = merged_cd_dens[node_id,:3]
        location = permute_map_coord_to_pdb(location,mapc,mapr,maps)
        All_location.append([location[k]+new_origin[k] for k in range(3)])
    return All_location
def Convert_LDPcoord_To_Reallocation(Coordinate_List, map_info_list):
    #1st filter out all the coordinates in this cluster
    mapc, mapr, maps, origin, nxstart, nystart, nzstart = map_info_list
    nstart = [nxstart,nystart,nzstart]
    nstart = permute_ns_coord_to_pdb(nstart,mapc,mapr,maps)
    new_origin = [origin[k]+nstart[k] for k in range(3)]
    All_location = []
    for node_id in range(len(Coordinate_List)):
        node_id = int(node_id)
        location =Coordinate_List[node_id]
        location = permute_map_coord_to_pdb(location,mapc,mapr,maps)
        All_location.append([location[k]+new_origin[k] for k in range(3)])
    return All_location

def calculate_merge_point_density(point,prob_array):
    before_merge_data = point.merged_data #[merge_to_id,x,y,z, density,merge_status]# 1st column is the id points to the merged points, which is the member id
    after_merge_data = point.merged_cd_dens
    Output_Prob = np.zeros(len(after_merge_data))
    count_isolate=0
    for k in range(len(before_merge_data)):
        merge_to_id,x,y,z,density,merge_status=before_merge_data[k]
        x,y,z = int(x),int(y),int(z)
        if merge_to_id==k and merge_status==-1:
            count_isolate+=1
            continue
        while merge_status==-1:
            merge_to_id = int(merge_to_id)
            merge_to_id,_,_,_,_,merge_status=before_merge_data[merge_to_id]
        final_id = int(merge_status)
        Output_Prob[final_id]+=prob_array[x,y,z]
    point.merge_prob = Output_Prob
    print("in total %d isolated grid points"%count_isolate)
    Output_Prob = np.zeros(len(after_merge_data))
    for k in range(len(after_merge_data)):
        x,y,z,_ = after_merge_data[k]
        x,y,z = int(x),int(y),int(z)
        Output_Prob[k]=prob_array[x,y,z]
    point.point_prob = Output_Prob
    return point


def build_LDP(input_mrc,sugar_density, sugar_Nact,origin_map_path,save_path,ext_name,params,map_info_list,relax_LDP=False):
    #relax_LDP only allow mean shift limited distance and limited length
    mean_shift_path = os.path.join(save_path, ext_name+'_mean_shift')
    sugar_point = Points(params, sugar_Nact)
    if relax_LDP:
        if not os.path.exists(mean_shift_path + '_cd_relax.txt') \
                or not os.path.exists(mean_shift_path + '_dens_rleax.txt'):
            input_mrc.general_mean_shift(sugar_density,sugar_point, mean_shift_path,constriant=True)
        else:
            input_mrc.load_general_mean_shift(sugar_density,sugar_point, mean_shift_path,constriant=True)
    else:
        if not os.path.exists(mean_shift_path + '_cd.txt') \
                or not os.path.exists(mean_shift_path + '_dens.txt'):
            input_mrc.general_mean_shift(sugar_density,sugar_point, mean_shift_path)
        else:
            input_mrc.load_general_mean_shift(sugar_density,sugar_point, mean_shift_path)
    #change the density to common value without dividing the Nori
    sugar_point.recover_density()
    sugar_point_path = os.path.join(save_path, ext_name+'_point.txt')
    # init_id, x,y,z,density, merged_to_id
    # can use the x,y,z here to assign detailed probability for each LDP points.
    if not os.path.exists(sugar_point_path):
        sugar_point.Merge_point(input_mrc, sugar_point_path)  # You will get a merged point file here.
    else:
        sugar_point.load_merge(input_mrc, sugar_point_path)
    merged_cd_dens = np.loadtxt(sugar_point_path[:-4] + 'onlymerged.txt')
    LDP_save_path = os.path.join(save_path,ext_name+"_LDP.mrc")
    save_LDP_map(LDP_save_path, merged_cd_dens, origin_map_path)
    mapc, mapr, maps, origin, nxstart, nystart, nzstart = map_info_list
    All_location = Extract_LDP_coord(merged_cd_dens,mapc, mapr, maps, origin, nxstart, nystart, nzstart)
    graph_path = os.path.join(save_path, ext_name+"_LDP.pdb")
    Show_Graph_Connect(All_location, [], graph_path)
    #Get each LDP's sum probability values from its neighbors
    sugar_point = calculate_merge_point_density(sugar_point,sugar_density)
    #plot in b factor
    ldp_prob_path = os.path.join(save_path, ext_name+"_LDPdens.cif")
    Show_Bfactor_cif(ext_name+"_dens",All_location,ldp_prob_path,sugar_point.merged_cd_dens[:,3])
    #turns out density showed very good results its correlation with real phosphate positions
    return sugar_point






