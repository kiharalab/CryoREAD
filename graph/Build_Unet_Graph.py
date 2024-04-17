import shutil
import mrcfile
from data_processing.map_utils import process_map_data
import numpy as np
from structure.MRC import MRC
from structure.Tree import Tree
from graph.LDP_ops import build_LDP,Build_Base_LDP,Build_Baseall_LDP,prepare_all_sugar_location
from graph.Graph_ops import  construct_graph
import os
from ops.os_operation import  mkdir
from graph.path_utils import collect_all_searched_path
from ops.fasta_utils import read_fasta
from graph.Assign_ops import Assign_Base_Main_Path_sugar,Assign_PhoLDP_Sugarpath,Assign_Base_Main_Path
from graph.visualize_utils import Visualize_Path_Base,Visualize_LDP_Path
import pickle
from graph.geo_utils import Match_Sugar_Base_Location
from graph.DP_ops import greedy_assign_PS
from graph.assemble_ops import  build_collision_table,solve_assignment
from graph.assignment_ext import Extend_Solve_Assignment_SP_support
from graph.reassign_ops import reassign_basedonfrag,merge_assign_geo_seq

from data_processing.format_pdb import format_pdb,remove_op3_pdb
from data_processing.Gen_MaskDRNA_map import Gen_MaskDRNA_map,Gen_MaskProtein_map
from graph.geo_structure_modeling import Build_Atomic_Structure

def refine_structure_global(init_pdb_path,format_pdb_path,root_save_path,frag_collect_dir,
                     chain_prob,origin_map_path,refined_pdb_path,params,DNA_label=False):
    try:
        #remove this to remove dependencies of pymol
        #os.system("pymol -cq ops/save_formated_pdb.py "+str(init_pdb_path)+" "+str(format_pdb_path))
        format_pdb(init_pdb_path,format_pdb_path)
        # 5.4.0 prepare the mask map
        mask_map_path = os.path.join(root_save_path,"mask_map.mrc")
        Gen_MaskDRNA_map(chain_prob,origin_map_path,mask_map_path,params['contour'],threshold=0.6)
        #run coot+phenix refinement
        output_dir = os.path.join(root_save_path,"Output")
        mkdir(output_dir)
        refine0_pdb_path = os.path.join(output_dir,"Refine_cycle0.pdb")
        shutil.copy(format_pdb_path,refine0_pdb_path)
        # 5.4 build final atomic structure with phenix.real_space_refine
        if params['colab']:
            os.system('cd %s; /content/phenix/phenix-1.20.1-4487/build/bin/phenix.real_space_refine %s %s '
                      'resolution=%.4f output.suffix="_phenix_refine" skip_map_model_overlap_check=True'%(frag_collect_dir,format_pdb_path,mask_map_path,params['resolution']))
        else:
            os.system('cd %s; phenix.real_space_refine %s %s resolution=%.4f output.suffix="_phenix_refine" skip_map_model_overlap_check=True'%(frag_collect_dir,format_pdb_path,mask_map_path,params['resolution']))
        gen_pdb_path = format_pdb_path[:-4]+"_phenix_refine_000.pdb"
        count_check=0
        while not os.path.exists(gen_pdb_path) and count_check<5:
            gen_pdb_path = format_pdb_path[:-4]+"_phenix_refine_00%d.pdb"%(count_check+1)
            count_check+=1
        if not os.path.exists(gen_pdb_path):
            print("please check final non-refined atomic structure in %s"%refine0_pdb_path)
            return
        try:
            #run coot+phenix refinement

            refine1_pdb_path = os.path.join(output_dir,"Refine_cycle1.pdb")
            remove_op3_pdb(gen_pdb_path,refine1_pdb_path)
            #do coot refinement
            refine2_pdb_path = os.path.join(output_dir,"Refine_cycle2.pdb")
            from coot.coot_refine_structure import coot_refine_structure
            if params['colab']:
                coot_software="/content/coot/bin/coot"
            else:
                coot_software="coot"
            coot_refine_structure(refine1_pdb_path,mask_map_path,refine2_pdb_path,coot_software)
            refine3_pdb_path = os.path.join(output_dir,"Refine_cycle3.pdb")
            if params['colab']:
                os.system('cd %s; /content/phenix/phenix-1.20.1-4487/build/bin/phenix.real_space_refine %s %s '
                          'resolution=%.4f output.suffix="_phenix_refine skip_map_model_overlap_check=True"'%(output_dir,refine2_pdb_path,mask_map_path,params['resolution']))
            else:
                os.system('cd %s; phenix.real_space_refine %s %s resolution=%.4f '
                          'output.suffix="_phenix_refine" skip_map_model_overlap_check=True'%(output_dir,refine2_pdb_path,mask_map_path,params['resolution']))
            phenix_final_pdb = refine2_pdb_path[:-4]+"_phenix_refine_000.pdb"
            count_check=0
            while not os.path.exists(phenix_final_pdb) and count_check<5:
                phenix_final_pdb = refine2_pdb_path[:-4]+"_phenix_refine_00%d.pdb"%(count_check+1)
                count_check+=1
            if os.path.exists(phenix_final_pdb):
                #shutil.move(phenix_final_pdb,refine3_pdb_path)
                format_pdb(phenix_final_pdb,refine3_pdb_path,DNA_label)
                print("please check final refined atomic structure in %s"%refine3_pdb_path)
                print("You can also check other refined output here %s"%output_dir)
            else:
                print("please check final refined atomic structure (pdb format) in this directory %s"%output_dir)
        except:
            #Final_Assemble_20_2_20_formated_phenix_refine_000.pdb
            print("refinement failed!")
            if os.path.exists(gen_pdb_path):
                #shutil.copy(gen_pdb_path,refined_pdb_path)
                format_pdb(gen_pdb_path,refined_pdb_path,DNA_label)
                print("please check final refined atomic structure in %s"%refined_pdb_path)
            else:
                print("please check final refined atomic structure in this directory %s"%frag_collect_dir)


    except Exception as e:
        print("Refinement failed because of the possible error",e)

def Build_Unet_Graph(origin_map_path,chain_prob_path,fasta_path,save_path,
                     gaussian_bandwidth,dcut,rdcut, params):
    root_save_path = os.path.split(save_path)[0]
    pho_prob_threshold=0.1#make sure all necessary edges are connected
    sugar_prob_threshold = 0.1
    base_prob_threshold = 0.25
    #0.1 read sequence information
    map_data, mapc, mapr, maps, origin, nxstart, nystart, nzstart = process_map_data(origin_map_path)
    map_info_list=[mapc, mapr, maps, origin, nxstart, nystart, nzstart]
    chain_class =8
    #["sugar", "phosphate","A","UT","C","G","protein","base"]
    #0.2 read sequence information
    if fasta_path is not None and os.path.exists(fasta_path) and os.path.getsize(fasta_path)>0:
        chain_dict,DNA_Label= read_fasta(input_fasta_path=fasta_path,dna_check=True)
        #DNA_Label = read_dna_label(chain_dict)
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
    sugar_graph,sugar_coordinate_list,sugar_edge_pairs,sugar_edge_d_dens = construct_graph(input_mrc,input_mrc.sugar_dens,sugar_point,chain_prob[0],chain_prob[-1],save_path,"sugar_graph",map_info_list,params,prob_threshold=sugar_prob_threshold,extend=True)

    #here edge density distance/probability since we prefer low distance+high prob
    #can't search once, it's too heavy loaded, hence divided into several graphs to calculate

    #2. Path searching with sugar graphs
    #2.1 divide subgraphs
    #mst = Tree(params)
    #sort edge, finding connections, label mst_label
    #connect_cid = mst.Setup_Connection(sugar_graph)
    #update with more relaxed subgraph construction, fully based on distance constraints
    #connect_cid = mst.Setup_Relaxed_Connection(sugar_graph)
    #subgraphs = sugar_graph.build_subgraph(connect_cid)
    #print("in total we have %d graphs"%len(subgraphs))
    #to allow more searches to avoid some failure
    subgraphs = [[k for k in range(sugar_graph.Nnode)]]

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
    save_seq_path = os.path.join(save_path,"Output_Structure_seq")
    save_path = os.path.join(save_path,"Output_Structure_noseq")
    mkdir(save_path)




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

    print("only apply geometry constraints for dynamic programming")
    greedy_save_path = os.path.join(save_path,"DP_geo_match")
    mkdir(greedy_save_path)
    from graph.DP_geo_ops import greedy_assign_geo
    overall_dict,frag_location_dict=greedy_assign_geo(All_Base_Path_List_sugar,All_Path_List,
    Path_P_align_list,Path_P_reverse_align_list,Pho_Prob_Refer_Dict,Pho_Prob_Refer_Reverse_Dict,
    sugar_point,map_info_list,greedy_save_path)

    #4.4 build collision table for assemble overlapped fragments
    frag_save_path = os.path.join(save_path,"AssembleFactory_geo")
    mkdir(frag_save_path)
    # 5 Final Atomic Structure Modeling

    # 5.1  prepare atom positions by paths
    sugar_path_based_location = prepare_all_sugar_location(All_Path_List,sugar_point,map_info_list)


    # 5.2 build initial atomic structure

    frag_collect_dir = os.path.join(frag_save_path,"atomic_geo")
    mkdir(frag_collect_dir)
    #from graph.geo_structure_modeling import Build_Atomic_Structure
    Build_Atomic_Structure(overall_dict,
    sugar_path_based_location, frag_collect_dir,
    Path_P_align_list,Path_P_reverse_align_list,pho_point,map_info_list,sugar_base_match_dict)
    init_pdb_path = os.path.join(frag_collect_dir,"Final_Assemble_geo.pdb")
    format_pdb_path = os.path.join(frag_collect_dir,"Final_Assemble_geo_formated.pdb")
    refined_pdb_path = os.path.join(save_path,"Final_Refinedgeo.pdb")
    if params['no_seqinfo'] or chain_dict is None:
        if not params['no_seqinfo'] and chain_dict is None:
             print("!!!parsing fasta input failed, we can not output structures considering sequence assignment!!!")
        # 5.3 reformat pdb for phenix to do refinement (including the last column in pdb file indicate atom type)
        if params['refine']:
            refine_structure_global(init_pdb_path,format_pdb_path,root_save_path,frag_collect_dir,
                     chain_prob,origin_map_path,refined_pdb_path,params,DNA_label=DNA_Label)
            final_pdb_path =  os.path.join(root_save_path,"CryoREAD.pdb")
            from ops.os_operation import collect_refine_pdb
            collect_refine_pdb(os.path.join(root_save_path,"Output"),refined_pdb_path,final_pdb_path)
        else:
            #output_dir = os.path.join(root_save_path,"Output")
            #mkdir(output_dir)
            nonrefined_pdb_path = os.path.join(root_save_path,"CryoREAD_norefine.pdb")
            #shutil.copy(init_pdb_path,nonrefined_pdb_path)
            format_pdb(init_pdb_path,nonrefined_pdb_path,DNA_Label)
            print("please check final output atomic structure in %s"%nonrefined_pdb_path)
        return
    else:
        format_pdb(init_pdb_path,format_pdb_path,DNA_Label)
        nonrefined_pdb_path = os.path.join(root_save_path,"CryoREAD_noseq.pdb")
        shutil.copy(format_pdb_path,nonrefined_pdb_path)
    overall_geo_dict= overall_dict
    frag_geo_location_dict = frag_location_dict
    save_path = save_seq_path
    mkdir(save_path)
    greedy_save_path = os.path.join(save_path,"DP_search_"+str(ldp_size)+"_"+str(checking_stride)+"_"+str(top_select))
    mkdir(greedy_save_path)
    # 6.1 greedy assign based on geomery as well as the
    if chain_dict is not None:
        if params['thread']==1:
            overall_dict,frag_location_dict=greedy_assign_PS(All_Base_Path_List_sugar,All_Path_List,
            Path_P_align_list,Path_P_reverse_align_list,Pho_Prob_Refer_Dict,Pho_Prob_Refer_Reverse_Dict,
            chain_prob,ldp_size,sugar_point,map_info_list,chain_dict,greedy_save_path,top_select,checking_stride)
        else:
            from graph.DP_ops import greedy_assign_PS_effective
            overall_dict,frag_location_dict=greedy_assign_PS_effective(All_Base_Path_List_sugar,All_Path_List,
            Path_P_align_list,Path_P_reverse_align_list,Pho_Prob_Refer_Dict,Pho_Prob_Refer_Reverse_Dict,
            chain_prob,ldp_size,sugar_point,map_info_list,chain_dict,greedy_save_path,top_select,checking_stride,
                                                                       num_cpus=params['thread'])
    else:
        print("!!!parsing fasta input failed, we can not output structures considering sequence assignment!!!")
        exit()

    # 6.2 Assemble fragments assignment together

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
    time_use = 3600*(len(sugar_point.merged_cd_dens)/1000)
    time_use = min(time_use,3600*10)
    if not os.path.exists(cur_final_assemble_path):

        solve_frag_combine_list = solve_assignment(collision_table,order_key_index,
                                                   order_chain_index,overall_dict,time_use=time_use)
        if len(solve_frag_combine_list)==0 and params['rule_soft']==1:
            print("no possible solution for assembling")
            print("please make contact with the developer for further help!")
            #return
        if len(solve_frag_combine_list)==0 and params['rule_soft']==0:
            print("no possible solution for assembling with hard rules")
            print("we will try to reassign for via soft rules")
            frag_save_path = os.path.join(save_path,"AssembleFactory_"+str(ldp_size)+"_"+str(checking_stride)+"_"+str(top_select))
            mkdir(frag_save_path)
            params['rule_soft']=1
            cur_final_assemble_path = os.path.join(frag_save_path,"assemble_frag_%d_%d_%d_soft%d.txt"%(ldp_size,checking_stride,top_select,params['rule_soft']))

            cur_collision_path = os.path.join(frag_save_path,"Collision_Table_%d_%d_%d_soft_%d.npy"%(ldp_size,checking_stride,top_select,params['rule_soft']))
            collision_table,order_key_index,order_chain_index,key_order_index=build_collision_table(All_Base_Path_List_sugar, checking_stride,
                        ldp_size,overall_dict,cur_collision_path,params['rule_soft'])
            solve_frag_combine_list = solve_assignment(collision_table,order_key_index,
                                                       order_chain_index,overall_dict,time_use=time_use)
        if len(solve_frag_combine_list)==0:
            print("no possible solution for assembling after trying strict rules and soft rules")
            print("please make contact with the developer for further help!")
            #return
        else:
            np.savetxt(cur_final_assemble_path,np.array(solve_frag_combine_list))
    else:
        solve_frag_combine_list = np.loadtxt(cur_final_assemble_path)
    print("loading solved possible fragment combination finished")

    #if we do not find any possible solution, we will directly call refine for the CryoREAD_noseq.pdb
    if len(solve_frag_combine_list)==0:
        print("*"*100)
        print("no solution find, use the noseq version as final output!")
        print("*"*100)
        final_pdb_path =  os.path.join(root_save_path,"CryoREAD.pdb")
        if not params['refine']:
            
            with open(final_pdb_path,"w") as f:
                f.write("#CryoREAD no_seq pdb, sequence information is not used because of unsolvable assembling!\n")
                with open(nonrefined_pdb_path,"r") as f2:
                    for line in f2:
                        f.write(line)
        else:
            from graph.refine_structure import refine_structure
            output_dir = os.path.join(root_save_path,"Output")
            os.makedirs(output_dir,exist_ok=True)
            refine_structure(nonrefined_pdb_path,origin_map_path,output_dir,params)
            from ops.os_operation import collect_refine_pdb
            refined_pdb_path =  os.path.join(root_save_path,"CryoREAD_tmp.pdb")
            collect_refine_pdb(output_dir,nonrefined_pdb_path,refined_pdb_path)
            with open(final_pdb_path,"w") as f:
                f.write("#CryoREAD no_seq pdb, sequence information is not used because of unsolvable assembling!\n")
                with open(refined_pdb_path,"r") as f2:
                    for line in f2:
                        f.write(line)
        
        return

    # 7.3 reassign for those overlapped regions
    frag_collect_dir = os.path.join(frag_save_path,"atomic_reassign_%d_%d_%d"%(ldp_size,checking_stride,top_select))
    mkdir(frag_collect_dir)
    overall_reassign_dict = reassign_basedonfrag(solve_frag_combine_list,order_key_index,order_chain_index,overall_dict,
        sugar_path_based_location, ldp_size,frag_collect_dir,checking_stride,top_select,chain_dict,
        All_Base_Path_List_sugar,All_Path_List,
        Path_P_align_list,Path_P_reverse_align_list,Pho_Prob_Refer_Dict,Pho_Prob_Refer_Reverse_Dict)

    # build a sequence-based only structure here
    from graph.structure_modeling import Build_Atomic_Model_nonoverlap_frag
    frag_collect_dir = os.path.join(frag_save_path,"atomic_seq")
    mkdir(frag_collect_dir)
    Extra_Added_Assign_Dict=Extend_Solve_Assignment_SP_support(frag_save_path,All_Path_List,solve_frag_combine_list,
        order_key_index,order_chain_index,overall_dict,define_ldp_size=ldp_size)
    print("We added %d assignment that have collision to fill those gaps"%(len(Extra_Added_Assign_Dict)))
    Build_Atomic_Model_nonoverlap_frag(overall_reassign_dict,
        sugar_path_based_location, ldp_size,frag_collect_dir,checking_stride,top_select,
        Path_P_align_list,Path_P_reverse_align_list,pho_point,map_info_list,
                    Extra_Added_Assign_Dict,sugar_base_match_dict,DNA_Label)
    init_pdb_path = os.path.join(frag_collect_dir,"Final_Assemble_%d_%d_%d.pdb"%(ldp_size,checking_stride,top_select))
    format_pdb_path = os.path.join(frag_collect_dir,"Final_Assemble_%d_%d_%d_formated.pdb"%(ldp_size,checking_stride,top_select))
    format_pdb(init_pdb_path,format_pdb_path,DNA_Label)
    output_dir = os.path.join(root_save_path,"Output")
    mkdir(output_dir)
    nonrefined_pdb_path = os.path.join(output_dir,"CryoREAD_seqonly.pdb")
    shutil.copy(format_pdb_path,nonrefined_pdb_path)

    # 7.4 merge dp assignment and the geo assignment
    final_assign_dict = merge_assign_geo_seq(overall_geo_dict,overall_reassign_dict)
    # 7.5 build atomic structure
    frag_collect_dir = os.path.join(frag_save_path,"atomic_all")
    mkdir(frag_collect_dir)

    Build_Atomic_Structure(final_assign_dict,
    sugar_path_based_location, frag_collect_dir,
    Path_P_align_list,Path_P_reverse_align_list,pho_point,map_info_list,sugar_base_match_dict)
    # 5.3 reformat pdb for phenix to do refinement (including the last column in pdb file indicate atom type)
    init_pdb_path = os.path.join(frag_collect_dir,"Final_Assemble_geo.pdb")
    format_pdb_path = os.path.join(frag_collect_dir,"Final_Assemble_seq_formated.pdb")
    refined_pdb_path = os.path.join(save_path,"Final_Refinedseq.pdb")
    if params['refine']:
        refine_structure_global(init_pdb_path,format_pdb_path,root_save_path,frag_collect_dir,
                 chain_prob,origin_map_path,refined_pdb_path,params,DNA_label=DNA_Label)
        final_pdb_path =  os.path.join(root_save_path,"CryoREAD.pdb")
        from ops.os_operation import collect_refine_pdb
        collect_refine_pdb(os.path.join(root_save_path,"Output"),refined_pdb_path,final_pdb_path)
    else:
        #output_dir = os.path.join(root_save_path,"Output")
        #mkdir(output_dir)
        nonrefined_pdb_path = os.path.join(root_save_path,"CryoREAD_norefine.pdb")
        #shutil.copy(init_pdb_path,nonrefined_pdb_path)
        format_pdb(init_pdb_path,nonrefined_pdb_path,DNA_Label)
        print("please check final output atomic structure in %s"%nonrefined_pdb_path)



