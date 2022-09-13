from graph.LDP_ops import permute_point_coord_to_global_coord,permute_global_coord_to_point_coord


_base_atoms = ["OP3","P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'"]
residue_atoms = {
    'A': ["N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4"],
    'G':  ["N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4"],
    'C': ["N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"],
    'U': ["N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"],
}
_pho_atoms=["OP3", "P", "OP1", "OP2", "O5'",]
_sugar_atoms = [ "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'"]
_backbone_atoms=['P',"O5'","C5'",  "C4'","C3'",  "O3'" ]
import numpy as np
def calculate_point_distance(a,b):

    dist= np.linalg.norm(np.array(a)-np.array(b))
    return dist
def calculate_cosine_angle(sugar_location,point_pho_location,point_next_pho_location):
    dist1 = calculate_point_distance(sugar_location,point_pho_location)
    dist2 = calculate_point_distance(point_pho_location,point_next_pho_location)
    dist3 = calculate_point_distance(sugar_location,point_next_pho_location)
    cos_value = (dist1**2+dist2**2-dist3**2)/(2*dist1*dist2)
    return cos_value
def read_refer_dict(pdb_path):
    pho_origin_dict = {}
    sugar_origin_dict = {}
    base_origin_dict = {}
    base_location_list=[]
    all_location_dict ={}
    with open(pdb_path,'r') as file:
        for line in file:
            if len(line)>4 and line[:4]=="ATOM":
                atom_name=line[13:16]
                atom_name = atom_name.replace(" ","")
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                location = [x,y,z]
                if atom_name in _pho_atoms:
                    pho_origin_dict[atom_name]=location
                elif atom_name in _sugar_atoms:
                    sugar_origin_dict[atom_name]=location
                else:
                    base_origin_dict[atom_name]=location
                    base_location_list.append(location)
                all_location_dict[atom_name]=location
    backbone_location_list =[]
    for atom_name in _backbone_atoms:
        if atom_name in pho_origin_dict:
            backbone_location_list.append(pho_origin_dict[atom_name])
        elif atom_name in sugar_origin_dict:
            backbone_location_list.append(sugar_origin_dict[atom_name])
    #calculate direction vector
    phosphate_location = np.array(pho_origin_dict['P'])

    backbone_location_list = np.array(backbone_location_list)
    base_center_location = np.array(base_location_list)
    base_center_location = np.mean(base_center_location,axis=0)
    bpp_cos_values = calculate_cosine_angle(base_center_location,backbone_location_list[0],backbone_location_list[-1])
    bp1_distance = calculate_point_distance(base_center_location,backbone_location_list[0])
    p1_p2_distance = calculate_point_distance(backbone_location_list[0],backbone_location_list[-1])
    p1_p2_vector = backbone_location_list[-1]-backbone_location_list[0]
    orthognal_location = backbone_location_list[0]+p1_p2_vector*(bpp_cos_values*bp1_distance)/p1_p2_distance
    p1p2_to_base_direction = base_center_location-orthognal_location
    orthognal_ratio = (bpp_cos_values*bp1_distance)/p1_p2_distance
    p1_to_base_direction = base_center_location- phosphate_location
    return all_location_dict, backbone_location_list,orthognal_ratio,p1p2_to_base_direction,p1_to_base_direction
standard_location_dict={}
standard_location_dict['A']=read_refer_dict("atomic/A.pdb")
standard_location_dict['C']=read_refer_dict("atomic/C.pdb")
standard_location_dict['G']=read_refer_dict("atomic/G.pdb")
standard_location_dict['U']=read_refer_dict("atomic/U.pdb")
standard_location_dict['DA']=read_refer_dict("atomic/DA.pdb")
standard_location_dict['DC']=read_refer_dict("atomic/DC.pdb")
standard_location_dict['DG']=read_refer_dict("atomic/DG.pdb")
standard_location_dict['DT']=read_refer_dict("atomic/DT.pdb")


def read_refer_dict_sugar(pdb_path):
    pho_origin_dict = {}
    sugar_origin_dict = {}
    base_origin_dict = {}
    base_location_list=[]
    all_location_dict ={}
    with open(pdb_path,'r') as file:
        for line in file:
            if len(line)>4 and line[:4]=="ATOM":
                atom_name=line[13:16]
                atom_name = atom_name.replace(" ","")
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                location = [x,y,z]
                if atom_name in _pho_atoms:
                    pho_origin_dict[atom_name]=location
                elif atom_name in _sugar_atoms:
                    sugar_origin_dict[atom_name]=location
                else:
                    base_origin_dict[atom_name]=location
                    base_location_list.append(location)
                all_location_dict[atom_name]=location
    backbone_location_list =[]
    for atom_name in _backbone_atoms:
        if atom_name in pho_origin_dict:
            backbone_location_list.append(pho_origin_dict[atom_name])
        elif atom_name in sugar_origin_dict:
            backbone_location_list.append(sugar_origin_dict[atom_name])
    backbone_location_list = np.array(backbone_location_list)
    sugar_location_list =[]
    for atom_name in sugar_origin_dict:
        sugar_location_list.append(sugar_origin_dict[atom_name])
    sugar_location_list = np.array(sugar_location_list)
    sugar_location_list = np.mean(sugar_location_list,axis=0)
    base_center_location = np.array(base_location_list)
    base_center_location = np.mean(base_center_location,axis=0)
    p1_p2_vector = backbone_location_list[-1]-backbone_location_list[0]
    next_sugar_location_list = sugar_location_list+p1_p2_vector


    bss_cos_values = calculate_cosine_angle(base_center_location,sugar_location_list,next_sugar_location_list)
    bs1_distance = calculate_point_distance(base_center_location,sugar_location_list)
    s1_s2_distance = calculate_point_distance(sugar_location_list,next_sugar_location_list)
    orthognal_location = sugar_location_list+p1_p2_vector*(bss_cos_values*bs1_distance)/s1_s2_distance
    s1s2_to_base_direction = base_center_location-orthognal_location
    orthognal_ratio = (bss_cos_values*bs1_distance)/s1_s2_distance
    s1_to_base_direction = base_center_location- sugar_location_list
    return all_location_dict, [sugar_location_list,next_sugar_location_list],orthognal_ratio,s1s2_to_base_direction,s1_to_base_direction


sugar_standard_location_dict={}
sugar_standard_location_dict['A']=read_refer_dict_sugar("atomic/A.pdb")
sugar_standard_location_dict['C']=read_refer_dict_sugar("atomic/C.pdb")
sugar_standard_location_dict['G']=read_refer_dict_sugar("atomic/G.pdb")
sugar_standard_location_dict['U']=read_refer_dict_sugar("atomic/U.pdb")
sugar_standard_location_dict['DA']=read_refer_dict_sugar("atomic/DA.pdb")
sugar_standard_location_dict['DC']=read_refer_dict_sugar("atomic/DC.pdb")
sugar_standard_location_dict['DG']=read_refer_dict_sugar("atomic/DG.pdb")
sugar_standard_location_dict['DT']=read_refer_dict_sugar("atomic/DT.pdb")



def rigid_transform_3D(A, B):
    assert len(A) == len(B)

    N = A.shape[0]  # total points
    centroid_A = np.mean(A, axis=0)
    centroid_B = np.mean(B, axis=0)

    # centre the points
    AA = A - np.tile(centroid_A, (N, 1))
    BB = B - np.tile(centroid_B, (N, 1))

    H = np.matmul(np.transpose(AA),BB)
    U, S, Vt = np.linalg.svd(H)
    R = np.matmul(Vt.T, U.T)

    # special reflection case
    if np.linalg.det(R) < 0:
        print("Reflection detected")
        Vt[2, :] *= -1
        R = np.matmul(Vt.T,U.T)

    t = -np.matmul(R, centroid_A) + centroid_B
    err = B - np.matmul(A,R.T) - t.reshape([1, 3])
    return R, t,err


def build_base_atomic_nuc_rigid(input_pho_position,next_pho_location, cur_base_location,nuc_type,map_info_list):
    point_pho_location = permute_global_coord_to_point_coord(input_pho_position,map_info_list)
    revert_back_pho_location = permute_point_coord_to_global_coord(point_pho_location,map_info_list)
    for k in range(3):
        assert abs(revert_back_pho_location[k]-input_pho_position[k])<=0.01
    point_next_pho_location = permute_global_coord_to_point_coord(next_pho_location,map_info_list)
    point_pho_location = np.array(point_pho_location)
    point_next_pho_location = np.array(point_next_pho_location)
    point_base_location = permute_global_coord_to_point_coord(cur_base_location,map_info_list)
    cos_base_p1_p2 = calculate_cosine_angle(point_base_location,point_pho_location,point_next_pho_location)
    base_p1_distance = calculate_point_distance(point_base_location,point_pho_location)
    p1_p2_distance = calculate_point_distance(point_pho_location,point_next_pho_location)
    orthognal_location = point_pho_location+(point_next_pho_location-point_pho_location)*(cos_base_p1_p2*base_p1_distance)/p1_p2_distance
    p1p2_base_orthognal_direction = point_base_location - orthognal_location

    current_standard_info = standard_location_dict[nuc_type]
    all_location_dict,standard_backbone_location_list, standard_orthognal_ratio,standard_p1p2_to_base_direction,standard_p1_base_direction = current_standard_info
    p1_p2_select_orthognal_location = point_pho_location + standard_orthognal_ratio*(point_next_pho_location-point_pho_location)
    base_direction_adjust_ratio = calculate_point_distance([0,0,0],standard_p1p2_to_base_direction)/calculate_point_distance([0,0,0],p1p2_base_orthognal_direction)
    p1_p2_adjust_base_location = p1_p2_select_orthognal_location+p1p2_base_orthognal_direction*base_direction_adjust_ratio

    standard_align_vector = []
    #standard_align_vector.append([0,0,0])
    #backbone coord
    #for k in range(len(standard_backbone_location_list)):
    #   standard_align_vector.append(standard_backbone_location_list[k]-standard_backbone_location_list[0])
    #add base coord

    #standard_align_vector.append(standard_backbone_location_list[-1]-standard_backbone_location_list[0])
    #standard_align_vector.append(standard_p1_base_direction)
    standard_align_vector.append(standard_backbone_location_list[0])
    standard_align_vector.append(standard_backbone_location_list[-1])
    standard_align_vector.append(standard_p1_base_direction+standard_backbone_location_list[0])
    for k in range(len(standard_align_vector)):
        standard_align_vector[k]=standard_align_vector[k]-standard_backbone_location_list[0]

    standard_align_vector = np.array(standard_align_vector)


    #build current vector
    current_bpp_vector = []#[[0,0,0]]
#   divide_number = len(standard_backbone_location_list)-1
#     p1_p2_direction = point_next_pho_location-point_pho_location
#     for k in range(1,divide_number+1):
#         current_bpp_vector.append(p1_p2_direction*k/divide_number)
#     #add base coord
#    # current_bpp_vector.append(point_next_pho_location-point_pho_location)
#     current_bpp_vector.append(p1_p2_adjust_base_location-point_pho_location)
    divide_number = len(standard_backbone_location_list)
    point_final_o_location = point_pho_location+(divide_number-1)/divide_number*(point_next_pho_location-point_pho_location)
    current_bpp_vector.append(point_pho_location)
    current_bpp_vector.append(point_final_o_location)
    current_bpp_vector.append(p1_p2_adjust_base_location)
    for k in range(len(current_bpp_vector)):
        current_bpp_vector[k]=current_bpp_vector[k]-point_pho_location
    current_bpp_vector = np.array(current_bpp_vector)

    rotation,translation,err = rigid_transform_3D(standard_align_vector,current_bpp_vector)
    all_locations = np.array(list(all_location_dict.values()))
    output_locations = np.matmul(all_locations, rotation.T) + translation.reshape([1, 3])
    search_key_list = list(all_location_dict.keys())
    pho_location_index = search_key_list.index("P")
    cur_pho_location = output_locations[pho_location_index]
    aligh_pho_shift = point_pho_location-cur_pho_location
    output_locations = output_locations+aligh_pho_shift
    output_locations = [permute_point_coord_to_global_coord(output_locations[k],map_info_list) for k in range(len(output_locations))]
    return search_key_list,output_locations
def build_base_atomic_nuc_rigid_sugar(input_pho_position,next_pho_location, cur_base_location,nuc_type,map_info_list):
    point_pho_location = permute_global_coord_to_point_coord(input_pho_position,map_info_list)
    revert_back_pho_location = permute_point_coord_to_global_coord(point_pho_location,map_info_list)
    for k in range(3):
        assert revert_back_pho_location[k]==input_pho_position[k]
    point_next_pho_location = permute_global_coord_to_point_coord(next_pho_location,map_info_list)
    point_pho_location = np.array(point_pho_location)
    point_next_pho_location = np.array(point_next_pho_location)
    point_base_location = permute_global_coord_to_point_coord(cur_base_location,map_info_list)
    cos_base_p1_p2 = calculate_cosine_angle(point_base_location,point_pho_location,point_next_pho_location)
    base_p1_distance = calculate_point_distance(point_base_location,point_pho_location)
    p1_p2_distance = calculate_point_distance(point_pho_location,point_next_pho_location)
    orthognal_location = point_pho_location+(point_next_pho_location-point_pho_location)*(cos_base_p1_p2*base_p1_distance)/p1_p2_distance
    #here actual is the sugar line to base direction
    p1p2_base_orthognal_direction = point_base_location - orthognal_location
    current_standard_info = sugar_standard_location_dict[nuc_type]
    all_location_dict,standard_backbone_location_list, standard_orthognal_ratio,standard_s1s2_to_base_direction,standard_s1_base_direction = current_standard_info

    s1_s2_select_orthognal_location = point_pho_location + standard_orthognal_ratio*(point_next_pho_location-point_pho_location)
    base_direction_adjust_ratio = calculate_point_distance([0,0,0],standard_s1s2_to_base_direction)/calculate_point_distance([0,0,0],p1p2_base_orthognal_direction)
    s1_s2_adjust_base_location = s1_s2_select_orthognal_location+p1p2_base_orthognal_direction*base_direction_adjust_ratio

    standard_align_vector = []
    standard_align_vector.append(standard_backbone_location_list[0])
    standard_align_vector.append(standard_backbone_location_list[-1])
    standard_align_vector.append(standard_s1_base_direction+standard_backbone_location_list[0])
    for k in range(len(standard_align_vector)):
        standard_align_vector[k]=standard_align_vector[k]-standard_backbone_location_list[0]

    standard_align_vector = np.array(standard_align_vector)


    #build current vector
    current_bpp_vector = []#[[0,0,0]]

    point_final_o_location = point_next_pho_location
    current_bpp_vector.append(point_pho_location)
    current_bpp_vector.append(point_final_o_location)
    current_bpp_vector.append(s1_s2_adjust_base_location)
    for k in range(len(current_bpp_vector)):
        current_bpp_vector[k]=current_bpp_vector[k]-point_pho_location
    current_bpp_vector = np.array(current_bpp_vector)

    rotation,translation,err = rigid_transform_3D(standard_align_vector,current_bpp_vector)
    all_locations = np.array(list(all_location_dict.values()))
    output_locations = np.matmul(all_locations, rotation.T) + translation.reshape([1, 3])
    search_key_list = list(all_location_dict.keys())
    pho_location_index = search_key_list.index("P")
    cur_pho_location = output_locations[pho_location_index]
    aligh_pho_shift = point_pho_location-cur_pho_location
    output_locations = output_locations+aligh_pho_shift
    output_locations = [permute_point_coord_to_global_coord(output_locations[k],map_info_list) for k in range(len(output_locations))]
    return search_key_list,output_locations
