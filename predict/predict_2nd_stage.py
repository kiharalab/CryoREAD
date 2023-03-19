from ops.os_operation import mkdir
import os
import mrcfile
import numpy as np
from data_processing.map_utils import permute_ns_coord_to_pdb,save_predict_specific_map,find_top_density
from predict.unet_detect_map_refine import unet_detect_map_refine
def Predict_2nd_Stage(input_map_path,prob_dir,model_path,save_path,
                      voxel_size,stride,batch_size,contour,params):
    #check if prediction exists, if it exists, then skip
    cur_predict_path = os.path.join(save_path,"Input")
    chain_predict_path = os.path.join(cur_predict_path,"chain_predictprob.npy")
    if os.path.exists(chain_predict_path) and os.path.getsize(chain_predict_path)>=1000:
        print("stage 2 prediction has been generated and saved in %s"%chain_predict_path)
        return

    #prob_dir: cascad 1st stage output directory
    cur_predict_path = os.path.join(prob_dir,"Input")
    chain_predict_path = os.path.join(cur_predict_path,"chain_predictprob.npy")
    base_predict_path = os.path.join(cur_predict_path,"base_predictprob.npy")
    base_predict_prob = np.load(base_predict_path)
    chain_predict_prob = np.load(chain_predict_path)
    #use map and contour information to filter background regions to save computation
    with mrcfile.open(input_map_path, permissive=True) as map_mrc:
         #normalize data
        map_data = np.array(map_mrc.data)
        # get the value serve as 1 in normalization
        map_data[map_data < 0] = 0
        print("map density range: %f %f"%(0,np.max(map_data)))
        percentile_98 = find_top_density(map_data,0.98)

        print("map hist log percentage 98: ",percentile_98)
        map_data[map_data > percentile_98] = percentile_98
        min_value = np.min(map_data)
        max_value = np.max(map_data)
        map_data = (map_data-min_value)/(max_value-min_value)
        nxstart, nystart, nzstart = map_mrc.header.nxstart, \
                                    map_mrc.header.nystart, \
                                    map_mrc.header.nzstart
        orig = map_mrc.header.origin
        orig = str(orig)
        orig = orig.replace("(", "")
        orig = orig.replace(")", "")
        orig = orig.split(",")
        nstart = [nxstart, nystart, nzstart]
        mapc = map_mrc.header.mapc
        mapr = map_mrc.header.mapr
        maps = map_mrc.header.maps
        print("detected mode mapc %d, mapr %d, maps %d" % (mapc, mapr, maps))
        nstart = permute_ns_coord_to_pdb(nstart, mapc, mapr, maps)
        new_origin = []
        for k in range(3):
            new_origin.append(float(orig[k]) + float(nstart[k]))

        print("Origin:", new_origin)
        train_save_path = os.path.join(save_path,"Input")
        mkdir(train_save_path)
        #adjust the contour level by the maximum value
        print("given contour %f"%contour)
        contour = contour/percentile_98
        print("revised contour %f"%contour)
        detection_all = unet_detect_map_refine(map_data,chain_predict_prob,base_predict_prob, model_path,
                                                voxel_size,stride,batch_size,
                                                    train_save_path,contour,params)

        chain_label_list = ["sugar", "phosphate","A","UT","C","G","protein","base"]
        for k,chain_name in enumerate(chain_label_list):
            cur_map_path = os.path.join(save_path, "chain_" + str(chain_name) + "_prob.mrc")
            save_predict_specific_map(cur_map_path, k , detection_all, input_map_path,label_only=False)





