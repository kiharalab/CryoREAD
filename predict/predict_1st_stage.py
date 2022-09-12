

from ops.os_operation import mkdir

import os
import mrcfile
import numpy as np
from data_processing.map_utils import save_predict_specific_map,permute_ns_coord_to_pdb,find_top_density
from predict.unet_detect_map_cascad import unet_detect_map_cascad


def Predict_1st_Stage(input_map_path,model_path,save_path,
                                 voxel_size,stride,batch_size,contour,params):
    #check if prediction exists, if it exists, then skip
    cur_predict_path = os.path.join(save_path,"Input")
    chain_predict_path = os.path.join(cur_predict_path,"chain_predictprob.npy")
    if os.path.exists(chain_predict_path) and os.path.getsize(chain_predict_path)>=1000:
        print("stage 1 prediction has been generated and saved in %s"%chain_predict_path)
        return
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
        detection_chain,detection_base = unet_detect_map_cascad(map_data,model_path,
                                                voxel_size,stride,batch_size,
                                                    train_save_path,contour,params)
        chain_label_list = ["sugar", "phosphate","base","protein",]
        for k,chain_name in enumerate(chain_label_list):
            cur_map_path = os.path.join(save_path, "chain_" + str(chain_name) + "_prob.mrc")
            save_predict_specific_map(cur_map_path, k , detection_chain, input_map_path,label_only=False)

        base_label_list = ["A", "UT","C","G",]
        for k,base_name in enumerate(base_label_list):
            cur_map_path = os.path.join(save_path, "base_" + str(base_name) + "_prob.mrc")
            save_predict_specific_map(cur_map_path, k , detection_base, input_map_path,label_only=False)

        #save extened base detection results
        base_region_detection = detection_chain[2]>0.5
        base_detect_compare = np.argmax(detection_base,axis=0)
        base_detect_compare[base_region_detection<=0]=-1

        for k,base_name in enumerate(base_label_list):
            #current_base= detection_base[k]>0.5
            current_win_base = base_detect_compare==k
            #combine_detection = current_base|current_win_base
            cur_map_path = os.path.join(save_path, "base_" + str(base_name) + "_win.mrc")
            save_predict_specific_map(cur_map_path, 1 ,current_win_base, input_map_path,label_only=True)

