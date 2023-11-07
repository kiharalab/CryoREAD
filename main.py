

import os
import mrcfile
import numpy as np
from ops.argparser import argparser
from ops.os_operation import mkdir
import time
def init_save_path(origin_map_path):
    save_path = os.path.join(os.getcwd(), 'Predict_Result')
    mkdir(save_path)
    map_name = os.path.split(origin_map_path)[1].replace(".mrc", "")
    map_name = map_name.replace(".map", "")
    map_name = map_name.replace("(","").replace(")","")
    save_path = os.path.join(save_path, map_name)
    mkdir(save_path)
    return save_path,map_name

if __name__ == "__main__":
    params = argparser()
    if params['mode']==0:
        #configure the running path
        running_dir = os.path.dirname(os.path.abspath(__file__))

        gpu_id = params['gpu']
        if gpu_id is not None:
            os.environ["CUDA_VISIBLE_DEVICES"] = gpu_id
        cur_map_path = os.path.abspath(params['F'])
        #process the map path if it's ending with .gz
        if ".gz"==cur_map_path[-3:]:
            from ops.os_operation import unzip_gz
            cur_map_path = unzip_gz(cur_map_path)
        model_dir = os.path.abspath(params['M'])
        if params['prediction_only'] or params['no_seqinfo']:
            fasta_path = None
        else:
            fasta_path = os.path.abspath(params['P'])
        if params['output'] is None:
            save_path,map_name = init_save_path(cur_map_path)
        else:
            save_path=params['output']
            map_name="input"
            mkdir(save_path)
        os.chdir(running_dir)
        from data_processing.Unify_Map import Unify_Map
        cur_map_path = Unify_Map(cur_map_path,os.path.join(save_path,map_name+"_unified.mrc"))
        from data_processing.Resize_Map import Resize_Map
        cur_map_path = Resize_Map(cur_map_path,os.path.join(save_path,map_name+".mrc"))
        #segment map to remove those useless regions
        from data_processing.map_utils import segment_map
        cur_new_map_path = os.path.join(save_path,map_name+"_segment.mrc")
        cur_map_path = segment_map(cur_map_path,cur_new_map_path,contour=0) #save the final prediction prob array space

        #do a pre check to notify user errors for contour
        with mrcfile.open(cur_map_path,permissive=True) as mrc:
            data=mrc.data
        if np.max(data)<=params['contour']:
            print("!!!Warning!!!Please provide contour level lower than maximum density value to run!")
            #exit()
            from data_processing.map_utils import automate_contour
            params['contour']=automate_contour(data)
            print("!!!Warning!!!reset contour level: %f"%params['contour'])

        from predict.predict_1st_stage import Predict_1st_Stage
        #1st stage cascad prediciton
        model_1st_stage_path = os.path.join(model_dir,"stage1_network.pth")
        save_path_1st_stage = os.path.join(save_path,"1st_stage_detection")
        mkdir(save_path_1st_stage)
        Predict_1st_Stage(cur_map_path,model_1st_stage_path,save_path_1st_stage,
                                 params['box_size'],params['stride'],params['batch_size'],params['contour'],params)
        #2nd stage refine prediction
        from predict.predict_2nd_stage import Predict_2nd_Stage
        model_2nd_stage_path = os.path.join(model_dir,"stage2_network.pth")
        save_path_2nd_stage = os.path.join(save_path,"2nd_stage_detection")
        mkdir(save_path_2nd_stage)
        Predict_2nd_Stage(cur_map_path,save_path_1st_stage,model_2nd_stage_path,save_path_2nd_stage,
                       params['box_size'],params['stride'],params['batch_size'],params['contour'],params)
        #0 gen protein map for other programs to run
        cur_predict_path = os.path.join(save_path_2nd_stage,"Input")
        chain_predict_path = os.path.join(cur_predict_path,"chain_predictprob.npy")
        chain_prob = np.load(chain_predict_path)#[sugar,phosphate,A,UT,C,G,protein,base]
        mask_map_path = os.path.join(save_path,"mask_protein.mrc")
        from data_processing.Gen_MaskDRNA_map import Gen_MaskProtein_map
        Gen_MaskProtein_map(chain_prob,cur_map_path,mask_map_path,params['contour'],threshold=0.3)
        if params['prediction_only']:
            print("Our prediction results are saved in %s with mrc format for visualization check."%save_path_2nd_stage)
            exit()
        #graph based atomic structure modeling
        gaussian_bandwidth = params['g'] #use 3
        dcut = params['m']#after meanshifting merge points distance<[float]
        rdcut = params['f']#remove ldp threshold 0.01
        from graph.Build_Unet_Graph import Build_Unet_Graph
        graph_save_path = os.path.join(save_path,"graph_atomic_modeling")
        mkdir(graph_save_path)


        Build_Unet_Graph(cur_map_path,chain_predict_path,fasta_path,graph_save_path,
                gaussian_bandwidth,dcut, rdcut,params)
    elif params['mode']==1:

        #structure refinement pipeline
        input_pdb = os.path.abspath(params['F'])
        input_map = os.path.abspath(params['M'])
        output_dir = os.path.abspath(params['P'])
        running_dir = os.path.dirname(os.path.abspath(__file__))
        os.chdir(running_dir)
        #resolution param also needed
        from graph.refine_structure import refine_structure
        refine_structure(input_pdb,input_map,output_dir,params)








