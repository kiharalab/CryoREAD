
from ops.os_operation import mkdir
import os
import mrcfile
import numpy as np

def permute_ns_coord_to_pdb(input_coord,mapc,mapr,maps):
    """
    :param input_coord: [x,y,z] coord from pdb
    :param mapc:
    :param mapr:
    :param maps:
    :return:
    """
    if mapc==1 and mapr==2 and maps==3:
        out_x = input_coord[0]#out_x coorespond to section
        out_y = input_coord[1]#out_y correspond to row
        out_z = input_coord[2]#out_z correspond to column
    elif mapc==1 and mapr==3 and maps==2:
        out_x = input_coord[0]
        out_y = input_coord[2]
        out_z = input_coord[1]
    elif mapc == 2 and mapr == 1 and maps == 3:
        out_x = input_coord[1]
        out_y = input_coord[0]
        out_z = input_coord[2]

    elif mapc == 2 and mapr == 3 and maps == 1:
        out_x = input_coord[2]
        out_y = input_coord[0]
        out_z = input_coord[1]
    elif mapc == 3 and mapr == 1 and maps == 2:
        out_x = input_coord[1]
        out_y = input_coord[2]
        out_z = input_coord[0]
    elif mapc == 3 and mapr == 2 and maps == 1:
        out_x = input_coord[2]
        out_y = input_coord[1]
        out_z = input_coord[0]
    else:
        exit()
    return [out_x, out_y, out_z]


def save_predict_specific_map(save_map_path,specific_class,prediction_array,
                              origin_map_path,label_only=False):
    prediction = np.array(prediction_array)
    #prediction[prediction!=specific_class]=0
    if label_only:
        prediction[prediction!=specific_class]=0
    else:
        prediction = prediction[specific_class]
    with mrcfile.open(origin_map_path, permissive=True) as mrc:
        prev_voxel_size = mrc.voxel_size
        prev_voxel_size_x = float(prev_voxel_size['x'])
        prev_voxel_size_y = float(prev_voxel_size['y'])
        prev_voxel_size_z = float(prev_voxel_size['z'])
        nx, ny, nz, nxs, nys, nzs, mx, my, mz = \
            mrc.header.nx, mrc.header.ny, mrc.header.nz, \
            mrc.header.nxstart, mrc.header.nystart, mrc.header.nzstart, \
            mrc.header.mx, mrc.header.my, mrc.header.mz
        orig = mrc.header.origin
        print("Origin:", orig)
        print("Previous voxel size:", prev_voxel_size)
        print("nx, ny, nz", nx, ny, nz)
        print("nxs,nys,nzs", nxs, nys, nzs)
        print("mx,my,mz", mx, my, mz)

        data_new = np.float32(prediction)
        mrc_new = mrcfile.new(save_map_path, data=data_new, overwrite=True)
        vsize = mrc_new.voxel_size
        vsize.flags.writeable = True
        vsize.x = 1.0
        vsize.y = 1.0
        vsize.z = 1.0
        mrc_new.voxel_size = vsize
        mrc_new.update_header_from_data()
        mrc_new.header.nxstart = nxs * prev_voxel_size_x
        mrc_new.header.nystart = nys * prev_voxel_size_y
        mrc_new.header.nzstart = nzs * prev_voxel_size_z
        mrc_new.header.mapc = mrc.header.mapc
        mrc_new.header.mapr = mrc.header.mapr
        mrc_new.header.maps = mrc.header.maps
        mrc_new.header.origin = orig
        mrc_new.update_header_stats()
        mrc.print_header()
        mrc_new.print_header()
        mrc_new.close()
        del data_new

#finding the maximum density used to normalize
def find_top_density(map_data,threshold):
    use_density=map_data[map_data>0]

    hist,bin_edges=np.histogram(use_density, bins=200)
    #print(hist)
    log_hist = [np.log(x) if x>0 else 0 for x in hist]
    sum_cutoff=np.sum(log_hist)*threshold
    cumulative=0
    for j in range(len(log_hist)):
        cumulative+=log_hist[j]
        if cumulative>sum_cutoff:
            return bin_edges[j]
    return bin_edges[-1]

def process_map_data(input_map_path):
    with mrcfile.open(input_map_path, permissive=True) as mrc:
        orig = mrc.header.origin
        orig = str(orig)
        orig = orig.replace("(","")
        orig = orig.replace(")", "")
        orig = orig.split(",")
        new_origin = []
        for k in range(3):
            new_origin.append(float(orig[k]))
        print("Origin:", new_origin)
        data = mrc.data
        mapc = mrc.header.mapc
        mapr = mrc.header.mapr
        maps = mrc.header.maps
        print("detected mode mapc %d, mapr %d, maps %d" % (mapc, mapr, maps))
        nxstart, nystart, nzstart = mrc.header.nxstart, mrc.header.nystart, mrc.header.nzstart
        #mapc to x, mapr to y maps to z
    return data, mapc, mapr, maps,new_origin,nxstart,nystart,nzstart


def permute_map_coord_to_pdb(input_coord,mapc,mapr,maps):
    """
    :param input_coord: [x,y,z] coord from pdb
    :param mapc:
    :param mapr:
    :param maps:
    :return:
    """
    if mapc==1 and mapr==2 and maps==3:
        out_x = input_coord[2]#out_x coorespond to section
        out_y = input_coord[1]#out_y correspond to row
        out_z = input_coord[0]#out_z correspond to column
    elif mapc==1 and mapr==3 and maps==2:
        out_x = input_coord[2]
        out_y = input_coord[0]
        out_z = input_coord[1]
    elif mapc == 2 and mapr == 1 and maps == 3:
        out_x = input_coord[1]
        out_y = input_coord[2]
        out_z = input_coord[0]

    elif mapc == 2 and mapr == 3 and maps == 1:
        out_x = input_coord[0]
        out_y = input_coord[2]
        out_z = input_coord[1]
    elif mapc == 3 and mapr == 1 and maps == 2:
        out_x = input_coord[1]
        out_y = input_coord[0]
        out_z = input_coord[2]
    elif mapc == 3 and mapr == 2 and maps == 1:
        out_x = input_coord[0]
        out_y = input_coord[1]
        out_z = input_coord[2]
    else:
        exit()
    return out_x, out_y, out_z

def permute_pdb_coord_to_map(input_coord,mapc,mapr,maps):
    """
    :param input_coord: [x,y,z] coord from pdb
    :param mapc:
    :param mapr:
    :param maps:
    :return:
    """
    if mapc==1 and mapr==2 and maps==3:
        out_x = input_coord[2]#out_x coorespond to section
        out_y = input_coord[1]#out_y correspond to row
        out_z = input_coord[0]#out_z correspond to column
    elif mapc==1 and mapr==3 and maps==2:
        out_x = input_coord[1]
        out_y = input_coord[2]
        out_z = input_coord[0]
    elif mapc == 2 and mapr == 1 and maps == 3:
        out_x = input_coord[2]
        out_y = input_coord[0]
        out_z = input_coord[1]

    elif mapc == 2 and mapr == 3 and maps == 1:
        out_x = input_coord[0]
        out_y = input_coord[2]
        out_z = input_coord[1]
    elif mapc == 3 and mapr == 1 and maps == 2:
        out_x = input_coord[1]
        out_y = input_coord[0]
        out_z = input_coord[2]
    elif mapc == 3 and mapr == 2 and maps == 1:
        out_x = input_coord[0]
        out_y = input_coord[1]
        out_z = input_coord[2]
    else:
        exit()
    return out_x, out_y, out_z
