
import imp
import numpy as np
import os
import datetime
import time
import torch
import torch.nn as nn
from ops.Logger import AverageMeter,ProgressMeter
from data_processing.DRNA_dataset import  Single_Dataset2
from model.Small_Unet_3Plus_DeepSup import Small_UNet_3Plus_DeepSup
def gen_input_data(map_data,chain_prob,base_prob,voxel_size,stride,contour,train_save_path):
    from progress.bar import Bar
    scan_x, scan_y, scan_z = map_data.shape
    chain_classes =len(chain_prob)
    base_classes = len(base_prob)
    count_voxel = 0
    count_iter=0
    Coord_Voxel = []
    bar = Bar('Preparing Input: ', max=int(np.ceil(scan_x/stride)*np.ceil(scan_y/stride)*np.ceil(scan_z/stride)))

    for x in range(0, scan_x, stride):
        x_end = min(x + voxel_size, scan_x)
        for y in range(0, scan_y, stride):
            y_end = min(y + voxel_size, scan_y)
            for z in range(0, scan_z, stride):
                bar.next()
                count_iter+=1
                z_end = min(z + voxel_size, scan_z)
                if x_end < scan_x:
                    x_start = x
                else:
                    x_start = x_end - voxel_size
                    if x_start<0:
                        x_start=0
                if y_end < scan_y:
                    y_start = y
                else:
                    y_start = y_end - voxel_size
                    if y_start<0:
                        y_start=0
                if z_end < scan_z:
                    z_start = z
                else:
                    z_start = z_end - voxel_size
                    if z_start<0:
                        z_start=0
                #already normalized
                segment_map_voxel = map_data[x_start:x_end, y_start:y_end, z_start:z_end]
                if contour<=0:
                    meaningful_density_count = len(np.argwhere(segment_map_voxel>0))
                    meaningful_density_ratio = meaningful_density_count/float(voxel_size**3)
                    if meaningful_density_ratio<=0.001:
                       # print("meaningful density ratio %f of current box, skip it!"%meaningful_density_ratio)
                        continue
                else:
                    meaningful_density_count = len(np.argwhere(segment_map_voxel > contour))
                    meaningful_density_ratio = meaningful_density_count / float(voxel_size ** 3)
                    if meaningful_density_ratio <= 0.001:
                       # print("no meaningful density ratio of current box, skip it!")
                        continue
                segment_input_voxel = np.zeros([chain_classes+base_classes, voxel_size,voxel_size,voxel_size])
                segment_input_voxel[:chain_classes,:x_end-x_start,:y_end-y_start,:z_end-z_start]=chain_prob[:,x_start:x_end, y_start:y_end, z_start:z_end]
                segment_input_voxel[chain_classes:,:x_end-x_start,:y_end-y_start,:z_end-z_start]=base_prob[:,x_start:x_end, y_start:y_end, z_start:z_end]
                #check values in segment input_voxel
                #different classes >=0.5 number should be bigger than 0.001
                count_meaningful = (segment_input_voxel>0.5).sum()
                meaningful_density_ratio = count_meaningful/ float(voxel_size ** 3)
                if meaningful_density_ratio <= 0.001:
                    #print("no meaningful predictions of current box in 1st stage, skip it!")
                    continue


                cur_path = os.path.join(train_save_path,"input_"+str(count_voxel)+".npy")
                np.save(cur_path,segment_input_voxel)
                Coord_Voxel.append([x_start,y_start,z_start])
                count_voxel+=1
                #print("2nd stage: %.2f percent scanning finished"%(count_iter/(scan_x*scan_y*scan_z/(stride**3))))
    bar.finish()
    Coord_Voxel = np.array(Coord_Voxel)
    coord_path = os.path.join(train_save_path,"Coord.npy")
    np.save(coord_path,Coord_Voxel)
    print("in 2nd stage, in total we have %d boxes"%len(Coord_Voxel))
    return Coord_Voxel

from model.Small_Unet_3Plus_DeepSup import Small_UNet_3Plus_DeepSup
import gc
def make_predictions(test_loader,model,Coord_Voxel,voxel_size,overall_shape,num_classes,run_type=0):
    avg_meters = {'data_time': AverageMeter('data_time'),
                  'train_time': AverageMeter('train_time')}
    progress = ProgressMeter(
        len(test_loader),
        [avg_meters['data_time'],
         avg_meters['train_time'],
         ],
        prefix="#Eval:")
    model.eval()
    end_time = time.time()
    scan_x, scan_y, scan_z = overall_shape
    Prediction_Matrix = np.zeros([num_classes,overall_shape[0],overall_shape[1],overall_shape[2]])
    #Count_Matrix = np.zeros(overall_shape)
    with torch.no_grad():
        for batch_idx, data in enumerate(test_loader):
            # input, atom_target, nuc_target = data
            input, cur_index = data
            #print(input.shape)
            avg_meters['data_time'].update(time.time() - end_time, input.size(0))
            cur_id = cur_index.detach().cpu().numpy()#test_loader.dataset.id_list[cur_index.detach().numpy()]
            input = input.cuda()
            outputs = model(input)
            if run_type==2:
                final_output = torch.softmax(torch.sigmoid(outputs[0]),dim=1).detach().cpu().numpy()
            elif run_type==1:
                final_output = torch.sigmoid(outputs[0]).detach().cpu().numpy()
            else:
                final_output = torch.softmax(outputs[0],dim=1).detach().cpu().numpy()
            avg_meters['train_time'].update(time.time() - end_time, input.size(0))
            progress.display(batch_idx)
            for k in range(len(cur_id)):
                tmp_index = cur_id[k]
                x_start, y_start, z_start = Coord_Voxel[int(tmp_index)]
                x_end,y_end,z_end = x_start+voxel_size,y_start+voxel_size,z_start+voxel_size
                if x_end < scan_x:
                    x_start = x_start
                else:
                    x_end = scan_x
                    x_start = x_end - voxel_size
                    if x_start<0:
                        x_start=0
                if y_end < scan_y:
                    y_start = y_start
                else:
                    y_end = scan_y
                    y_start = y_end - voxel_size
                    if y_start<0:
                        y_start=0
                if z_end < scan_z:
                    z_start = z_start
                else:
                    z_end=scan_z
                    z_start = z_end - voxel_size
                    if z_start<0:
                        z_start=0
                #print(final_output[k].shape)
                #print(Prediction_Matrix[:,x_start:x_end,y_start:y_end,z_start:z_end].shape)
                #pred_label=np.argmax(final_output[k],axis=1)
                #count_positive= len(np.argwhere(pred_label!=0))
                #print("%d example with %d positive predictions"%(k,count_positive))
                Prediction_Matrix[:,x_start:x_end,y_start:y_end,z_start:z_end] =np.maximum(Prediction_Matrix[:,x_start:x_end,y_start:y_end,z_start:z_end],final_output[k][:,:x_end-x_start,:y_end-y_start,:z_end-z_start])

                #Count_Matrix[x_start:x_end,y_start:y_end,z_start:z_end]+=1
            if batch_idx%1000==0:
                for j in range(num_classes):
                    count_positive = len(np.argwhere(Prediction_Matrix[j]>=0.5))
                    print("%d classes already detected %d voxels"%(j,count_positive))
            end_time = time.time()
            del final_output
            #del pred_label
            del outputs
            del input
            gc.collect()
    #Prediction_Matrix = Prediction_Matrix/Count_Matrix
    #replace nan with 0
    Prediction_Matrix[np.isnan(Prediction_Matrix)] = 0
    Prediction_Label = np.argmax(Prediction_Matrix,axis=0)


    return Prediction_Matrix,Prediction_Label
def unet_detect_map_refine(map_data,chain_prob,base_prob,resume_model_path,voxel_size,
                    stride,batch_size,train_save_path,contour,params):
    coord_path = os.path.join(train_save_path, "Coord.npy")
    if os.path.exists(coord_path):
        Coord_Voxel = np.load(coord_path)
    else:
        Coord_Voxel = gen_input_data(map_data,chain_prob,base_prob,voxel_size, stride, contour,train_save_path)
    overall_shape = map_data.shape
    test_dataset = Single_Dataset2(train_save_path,"input_")
    test_loader = torch.utils.data.DataLoader(
        test_dataset,
        batch_size=batch_size,
        pin_memory=True,
        shuffle=False,
        num_workers=params['num_workers'],
        drop_last=False)
    chain_class = len(chain_prob)
    base_class = len(base_prob)
    output_classes = chain_class+base_class
    model = Small_UNet_3Plus_DeepSup(in_channels=chain_class+base_class,
                                     n_classes=output_classes,
                                     feature_scale=4,
                                     is_deconv=True,
                                     is_batchnorm=True)
    model = model.cuda()
    model = nn.DataParallel(model, device_ids=None)
    state_dict = torch.load(resume_model_path)
    msg = model.load_state_dict(state_dict['state_dict'])
    print("model loading: ",msg)
    cur_prob_path = os.path.join(train_save_path, "chain_predictprob.npy")
    cur_label_path = os.path.join(train_save_path, "chain_predict.npy")
    if os.path.exists(cur_prob_path) and os.path.exists(cur_label_path):
        Prediction_Matrix =np.load(cur_prob_path)
        Prediction_Label = np.load(cur_label_path)

    else:
        Prediction_Matrix,Prediction_Label = make_predictions(test_loader,model,Coord_Voxel,
                                            voxel_size,overall_shape,
                                            output_classes,run_type=1)#must sigmoid activated

        np.save(cur_prob_path,Prediction_Matrix)
        np.save(cur_label_path, Prediction_Label)
    return Prediction_Matrix





