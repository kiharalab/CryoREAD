
import numpy as np
import os
import datetime
import time
import torch
import torch.nn as nn
from ops.Logger import AverageMeter,ProgressMeter
from data_processing.DRNA_dataset import Single_Dataset
from model.Cascade_Unet import Cascade_Unet


def gen_input_data(map_data,voxel_size,stride,contour,train_save_path):
    scan_x, scan_y, scan_z = map_data.shape
    count_voxel = 0
    count_iter=0
    Coord_Voxel = []
    from progress.bar import Bar
    bar = Bar('Preparing Input: ', max=int(np.ceil(scan_x/stride)*np.ceil(scan_y/stride)*np.ceil(scan_z/stride)))


    for x in range(0, scan_x, stride):
        x_end = min(x + voxel_size, scan_x)
        for y in range(0, scan_y, stride):
            y_end = min(y + voxel_size, scan_y)
            for z in range(0, scan_z, stride):
                count_iter+=1
                bar.next()
                #print("1st stage: %.4f percent scanning finished"%(count_iter*100/(scan_x*scan_y*scan_z/(stride**3))),"location %d %d %d"%(x,y,z))
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
                segment_map_voxel = np.zeros([voxel_size,voxel_size,voxel_size])
                segment_map_voxel[:x_end-x_start,:y_end-y_start,:z_end-z_start]=map_data[x_start:x_end, y_start:y_end, z_start:z_end]
                if contour<=0:
                    meaningful_density_count = len(np.argwhere(segment_map_voxel>0))
                    meaningful_density_ratio = meaningful_density_count/float(voxel_size**3)
                    if meaningful_density_ratio<=0.001:
                        #print("no meaningful density ratio %f in current scanned box, skip it!"%meaningful_density_ratio)
                        continue
                else:
                    meaningful_density_count = len(np.argwhere(segment_map_voxel > contour))
                    meaningful_density_ratio = meaningful_density_count / float(voxel_size ** 3)
                    if meaningful_density_ratio <= 0.001:
                       # print("no meaningful density ratio in current scanned box, skip it!")
                        continue
                cur_path = os.path.join(train_save_path,"input_"+str(count_voxel)+".npy")
                np.save(cur_path,segment_map_voxel)
                Coord_Voxel.append([x_start,y_start,z_start])
                count_voxel+=1
    bar.finish()
    Coord_Voxel = np.array(Coord_Voxel)
    coord_path = os.path.join(train_save_path,"Coord.npy")
    np.save(coord_path,Coord_Voxel)
    print("In total we prepared %d boxes as input"%(len(Coord_Voxel)))
    return Coord_Voxel


def make_predictions(test_loader,model,Coord_Voxel,voxel_size,overall_shape,num_classes,base_classes):
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
    Base_Matrix = np.zeros([base_classes,overall_shape[0],overall_shape[1],overall_shape[2]])
    Count_Matrix = np.zeros(overall_shape)
    #average for overlap regions
    with torch.no_grad():
        for batch_idx, data in enumerate(test_loader):
            # input, atom_target, nuc_target = data
            input, cur_index = data
            #print(input.shape)
            avg_meters['data_time'].update(time.time() - end_time, input.size(0))
            cur_id = cur_index.detach().cpu().numpy()#test_loader.dataset.id_list[cur_index.detach().numpy()]
            input = input.cuda()
            chain_outputs, base_outputs  = model(input)
            final_output = torch.sigmoid(chain_outputs[0]).detach().cpu().numpy()
            final_base = torch.sigmoid(base_outputs[0]).detach().cpu().numpy()

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

                Prediction_Matrix[:,x_start:x_end,y_start:y_end,z_start:z_end] += final_output[k][:,:x_end-x_start,:y_end-y_start,:z_end-z_start]
                Base_Matrix[:,x_start:x_end,y_start:y_end,z_start:z_end] += final_base[k][:,:x_end-x_start,:y_end-y_start,:z_end-z_start]
                Count_Matrix[x_start:x_end,y_start:y_end,z_start:z_end]+=1
            if batch_idx%1000==0:
                for j in range(num_classes):
                    count_positive = len(np.argwhere(Prediction_Matrix[j]>=0.5))
                    print("%d classes already detected %d voxels"%(j,count_positive))
            end_time = time.time()
    Prediction_Matrix = Prediction_Matrix/Count_Matrix
    #replace nan with 0
    Prediction_Matrix[np.isnan(Prediction_Matrix)] = 0
    Prediction_Label = np.argmax(Prediction_Matrix,axis=0)

    Base_Matrix = Base_Matrix/Count_Matrix
    Base_Matrix[np.isnan(Base_Matrix)] = 0
    #New_Base_Matrix = np.zeros([base_classes-1,overall_shape[0],overall_shape[1],overall_shape[2]])
    New_Base_Matrix = Base_Matrix[:-1]
    New_Base_Matrix[1] = np.maximum(Base_Matrix[1],Base_Matrix[-1])#merge u and t predictions

    Base_Label = np.argmax(New_Base_Matrix,axis=0)
    return Prediction_Matrix,Prediction_Label,New_Base_Matrix,Base_Label



def unet_detect_map_cascad(map_data,resume_model_path,voxel_size,
                    stride,batch_size,train_save_path,contour,params):

    coord_path = os.path.join(train_save_path, "Coord.npy")
    if os.path.exists(coord_path):
        Coord_Voxel = np.load(coord_path)
    else:
        Coord_Voxel = gen_input_data(map_data,voxel_size, stride, contour,train_save_path)
    overall_shape = map_data.shape
    test_dataset = Single_Dataset(train_save_path,"input_")
    test_loader = torch.utils.data.DataLoader(
        test_dataset,
        batch_size=batch_size,
        shuffle=False,
        pin_memory=True,
        num_workers=params['num_workers'],
        drop_last=False)
    chain_class = 4
    base_class = 5
    model =  Cascade_Unet(in_channels=1,#include density and probability array
                           n_classes1=chain_class,
                          n_classes2= base_class,
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
    cur_baseprob_path = os.path.join(train_save_path, "base_predictprob.npy")
    cur_baselabel_path = os.path.join(train_save_path, "base_predict.npy")
    if os.path.exists(cur_prob_path) and os.path.exists(cur_label_path):
        Prediction_Matrix =np.load(cur_prob_path)
        Prediction_Label = np.load(cur_label_path)
        Base_Matrix = np.load(cur_baseprob_path)
        Base_Label  = np.load(cur_baselabel_path)
    else:
        Prediction_Matrix,Prediction_Label,Base_Matrix,Base_Label = make_predictions(test_loader,model,Coord_Voxel,
                                            voxel_size,overall_shape,
                                            chain_class,base_class)

        np.save(cur_prob_path,Prediction_Matrix)
        np.save(cur_label_path, Prediction_Label)
        np.save(cur_baseprob_path,Base_Matrix)
        np.save(cur_baselabel_path,Base_Label)
    #save disk space for the generated input boxes
    os.system("rm "+train_save_path+"/input*")
    return Prediction_Matrix,Base_Matrix

