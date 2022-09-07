import shutil

from data_processing.map_utils import process_map_data
import numpy as np
from structure.MRC import MRC
from structure.Tree import Tree
from graph.LDP_ops import build_LDP,Build_Base_LDP,Build_Baseall_LDP,prepare_all_sugar_location
from graph.Graph_ops import  construct_graph
import os
from ops.os_operation import  mkdir
from graph.path_utils import collect_all_searched_path
from ops.fasta_utils import read_fasta,read_dna_label
from graph.Assign_ops import Assign_Base_Main_Path_sugar,Assign_PhoLDP_Sugarpath,Assign_Base_Main_Path
from graph.visualize_utils import Visualize_Path_Base,Visualize_LDP_Path
import pickle
from graph.geo_utils import Match_Sugar_Base_Location
from graph.DP_ops import greedy_assign_PS
from graph.assemble_ops import  build_collision_table,solve_assignment
from graph.assignment_ext import Extend_Solve_Assignment_SP_support

from graph.structure_modeling import  Build_Atomic_Structure

def Build_Unet_Graph(origin_map_path,chain_prob_path,fasta_path,save_path,
                     gaussian_bandwidth,dcut,rdcut, params):
    pho_prob_threshold=0.1#make sure all necessary edges are connected
    sugar_prob_threshold = 0.1
    base_prob_threshold = 0.25
    #0.1 read sequence information
    map_data, mapc, mapr, maps, origin, nxstart, nystart, nzstart = process_map_data(origin_map_path)
    map_info_list=[mapc, mapr, maps, origin, nxstart, nystart, nzstart]
    chain_class =8
    #["sugar", "phosphate","A","UT","C","G","protein","base"]
    #0.2 read sequence information
    if os.path.exists(fasta_path) and os.path.getsize(fasta_path)>0:
        chain_dict = read_fasta(input_fasta_path=fasta_path)
        DNA_Label = read_dna_label(chain_dict)
        print("we have %d chains in provided fasta files"%(len(chain_dict)))
    else:
        chain_dict = None
        DNA_Label=False#default processing as RNA


    chain_prob = np.load(chain_prob_path)#[sugar,phosphate,A,UT,C,G,protein,base]
    input_mrc = MRC(origin_map_path, gaussian_bandwidth)

    #1. chain tracing
    sp_prob = chain_prob[0]+chain_prob[1]
    pho_prob = chain_prob[1]
    sugar_prob = chain_prob[0]

    input_mrc.upsampling_pho_prob(pho_prob,threshold=pho_prob_threshold,filter_array=None)
    input_mrc.upsampling_sugar_prob(sugar_prob,threshold=sugar_prob_threshold,filter_array=None)

    #1.1 LDP construction based on probability map
    pho_point_path = os.path.join(save_path,"pho_LDP")
    mkdir(pho_point_path)
    pho_point= build_LDP(input_mrc,input_mrc.pho_dens, input_mrc.pho_Nact,origin_map_path,pho_point_path,"pho",params,map_info_list)
    sugar_point_path = os.path.join(save_path,"sugar_LDP")
    mkdir(sugar_point_path)
    sugar_point = build_LDP(input_mrc,input_mrc.sugar_dens,input_mrc.sugar_Nact,origin_map_path,sugar_point_path,"sugar",params,map_info_list)

    #1.2 graph construction: edge constructions
    #pho_graph,pho_coordinate_list,pho_edge_pairs,pho_edge_d_dens = construct_graph(input_mrc,input_mrc.pho_dens, pho_point,chain_prob[1],chain_prob[-1],save_path,"pho_graph",map_info_list,params,prob_threshold=pho_prob_threshold,extend=True)
    sugar_graph,sugar_coordinate_list,sugar_edge_pairs,sugar_edge_d_dens = construct_graph(input_mrc,input_mrc.sugar_dens,sugar_point,chain_prob[0],chain_prob[-1],save_path,"sugar_graph",map_info_list,params,prob_threshold=sugar_prob_threshold,extend=False)

    #here edge density distance/probability since we prefer low distance+high prob
    #can't search once, it's too heavy loaded, hence divided into several graphs to calculate

    #2. Path searching with sugar graphs
    #2.1 divide subgraphs
    mst = Tree(params)
    #sort edge, finding connections, label mst_label
    connect_cid = mst.Setup_Connection(sugar_graph)
    subgraphs = sugar_graph.build_subgraph(connect_cid)
    print("in total we have %d graphs"%len(subgraphs))

    #2.2 collect all possible paths
    sugar_search_path = os.path.join(save_path,"sugar_search")
    mkdir(sugar_search_path)
    All_Path_List = collect_all_searched_path(subgraphs, sugar_point,sugar_graph,sugar_search_path,map_info_list,params)
    print("in total we collected %d path candidates"%len(All_Path_List))

    #3 Base Assignment
    #3.1 Cluster representative points of base LDPs
    #['A','UT','C','G','base'] ldps in the list
    base_save_path = os.path.join(save_path,"base_LDP")
    Base_LDP_List=Build_Base_LDP(input_mrc,chain_prob,base_prob_threshold,origin_map_path,params,map_info_list,base_save_path,filter_type=0)
    base_point= Build_Baseall_LDP(input_mrc,chain_prob,base_prob_threshold,origin_map_path,params,map_info_list,base_save_path)
    Base_LDP_List.append(base_point)

    #3.2 collect all possible paths
    All_Base_Path_List_sugar = Assign_Base_Main_Path_sugar(All_Path_List,sugar_point.merged_cd_dens[:,:3],
                    Base_LDP_List,chain_prob[2:6],cut_off_length=10)
    Visualize_Path_Base(os.path.join(save_path,"base_assignsugar_combine"),All_Path_List,
                        All_Base_Path_List_sugar,sugar_point.merged_cd_dens,
                        map_info_list,DNA_Label,sugar_visual=True)

    #3.3 define the assign correlation for sugar point and base ldp point.
    import pickle
    base_refer_path = os.path.join(save_path,"base_location_refer_dict.pkl")
    if not os.path.exists(base_refer_path):
        sugar_base_match_dict = Match_Sugar_Base_Location(All_Path_List,sugar_point,Base_LDP_List,map_info_list)
        with open(base_refer_path, 'wb') as handle:
            pickle.dump(sugar_base_match_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        with open(base_refer_path, 'rb') as handle:
            sugar_base_match_dict= pickle.load(handle)
    #4 DP assignment for each possible
    #4.1 DP assign for each fragments in the searched paths by applying dp for known sequences information.
    top_select = params['top_select']
    #mean value should be mean of all maps in reasonable regions.
    frag_size = params['frag_size']
    ldp_size = frag_size#int(frag_size*1.5)
    checking_stride = params['frag_stride']

    save_path = os.path.join(save_path,"Output_Structure")
    mkdir(save_path)

    greedy_save_path = os.path.join(save_path,"DP_search_"+str(ldp_size)+"_"+str(checking_stride)+"_"+str(top_select))
    mkdir(greedy_save_path)
    #4.1 align each phosphate to sugar to make combined DP
    #assign pho ldp locations to the path based on a projection to the sugar path, know which its previous node, which its next node
    #for each sugar node, assign its prev pointer to pho and next pointer to another pho
    Path_P_align_list,Path_P_reverse_align_list,Path_Pho_ID_List= Assign_PhoLDP_Sugarpath(All_Path_List,sugar_point,pho_point)
    #visualize the pho pdbs
    for k,path in enumerate(Path_Pho_ID_List):
        Visualize_LDP_Path(os.path.join(save_path,"support_p_%d.cif"%k),"support_p_%d"%k,path,pho_point.merged_cd_dens[:,:3],map_info_list)

    # 4.2 match pho-based base prob to each pho positions
    Pho_Prob_Refer_List,Pho_Prob_Refer_Dict = Assign_Base_Main_Path(Path_Pho_ID_List,pho_point.merged_cd_dens[:,:3], Base_LDP_List,chain_prob[2:6],cut_off_length=10,return_dict=True)

    Pho_Prob_Refer_List_Reverse,Pho_Prob_Refer_Reverse_Dict = Assign_Base_Main_Path(Path_Pho_ID_List,pho_point.merged_cd_dens[:,:3], Base_LDP_List,chain_prob[2:6],cut_off_length=10,reverse_flag=True,return_dict=True)

    #4.3 DP assign for each fragments in the searched paths by applying dp for known sequences information.



    if chain_dict is not None:
        overall_dict,frag_location_dict=greedy_assign_PS(All_Base_Path_List_sugar,All_Path_List,
        Path_P_align_list,Path_P_reverse_align_list,Pho_Prob_Refer_Dict,Pho_Prob_Refer_Reverse_Dict,
        chain_prob,ldp_size,sugar_point,map_info_list,chain_dict,greedy_save_path,top_select,checking_stride)


    else:
        print("only apply geometry constraints for dynamic programming")

    #4.4 build collision table for assemble overlapped fragments
    if params['rule_soft']==1:
        frag_save_path = os.path.join(save_path,"AssembleFactory_"+str(ldp_size)+"_"+str(checking_stride)+"_"+str(top_select))
    else:
        #pose strict rules for assembling
        frag_save_path = os.path.join(save_path,"AssembleFactory_strict_"+str(ldp_size)+"_"+str(checking_stride)+"_"+str(top_select))
    mkdir(frag_save_path)
    cur_final_assemble_path = os.path.join(frag_save_path,"assemble_frag_%d_%d_%d_soft%d.txt"%(ldp_size,checking_stride,top_select,params['rule_soft']))

    cur_collision_path = os.path.join(frag_save_path,"Collision_Table_%d_%d_%d_soft_%d.npy"%(ldp_size,checking_stride,top_select,params['rule_soft']))
    collision_table,order_key_index,order_chain_index,key_order_index=build_collision_table(All_Base_Path_List_sugar, checking_stride,
                        ldp_size,overall_dict,cur_collision_path,params['rule_soft'])

    #4.5 build assemble from pre-defined collision table

    if not os.path.exists(cur_final_assemble_path):

        solve_frag_combine_list = solve_assignment(collision_table,order_key_index,order_chain_index,overall_dict)
        np.savetxt(cur_final_assemble_path,np.array(solve_frag_combine_list))
    else:
        solve_frag_combine_list = np.loadtxt(cur_final_assemble_path)
    print("loading solved possible fragment combination finished")


    # 4.6 extend fragment assignment for regions without acceptable assignment
    Extra_Added_Assign_Dict=Extend_Solve_Assignment_SP_support(frag_save_path,All_Path_List,solve_frag_combine_list,
    order_key_index,order_chain_index,overall_dict,define_ldp_size=ldp_size)
    print("We added %d assignment that have collision to fill those gaps"%(len(Extra_Added_Assign_Dict)))


    # 5 Final Atomic Structure Modeling

    # 5.1  prepare atom positions by paths
    sugar_path_based_location = prepare_all_sugar_location(All_Path_List,sugar_point,map_info_list)


    # 5.2 build initial atomic structure
    frag_collect_dir = os.path.join(frag_save_path,"atomic_assemble_frag_%d_%d_%d"%(ldp_size,checking_stride,top_select))
    mkdir(frag_collect_dir)

    Build_Atomic_Structure(solve_frag_combine_list,order_key_index,order_chain_index,overall_dict,
    sugar_path_based_location, ldp_size,frag_collect_dir,checking_stride,top_select,
    Path_P_align_list,Path_P_reverse_align_list,pho_point,map_info_list,Extra_Added_Assign_Dict,
                           sugar_base_match_dict,DNA_Label)


    # 5.3 build final atomic structure with phenix.real_space_refine
    init_pdb_path = os.path.join(frag_save_path,"Final_Assemble_%d_%d_%d.pdb"%(ldp_size,checking_stride,top_select))
    refined_pdb_path = os.path.join(frag_save_path,"Final_Refined_%d_%d_%d.pdb"%(ldp_size,checking_stride,top_select))
    os.system('phenix.real_space_refine %s %s resolution=%.4f output.suffix="_phenix_refine"'%(init_pdb_path,origin_map_path,params['resolution']))
    gen_pdb_path = os.path.join(frag_save_path,"Final_Assemble_%d_%d_%dphenix_refine.pdb"%(ldp_size,checking_stride,top_select))
    if os.path.exists(gen_pdb_path):
        shutil.copy(gen_pdb_path,refined_pdb_path)
    print("please check final refined atomic structure in %s"%refined_pdb_path)



