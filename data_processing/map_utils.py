
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

def save_dens_map(save_map_path,new_dens,
                              origin_map_path):

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

        data_new = np.float32(new_dens)
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
def increase_map_density(input_path,output_path,add_contour):
    add_contour=abs(add_contour)
    with mrcfile.open(input_path,permissive=True) as mrc:
        data=mrc.data
        data=np.float32(data)
    data = data+add_contour
    save_dens_map(output_path,data,input_path)
    return output_path
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

def automate_contour(map_data):
    use_density=map_data[map_data>1e-6]
    #hist,bin_edges=np.histogram(use_density, bins=1000)
    sorted_array = np.sort(use_density)
    select_index=int(len(sorted_array)/1000) #usually 0
    return sorted_array[select_index]
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

def segment_map(input_map,output_map,contour=0):
    """
    segment meaningful region of a map
    :param input_map:
    :param output_map:
    :return:
    generate a new small size map
    """
    with mrcfile.open(input_map, permissive=True) as mrc:
        prev_voxel_size = mrc.voxel_size
        prev_voxel_size_x = float(prev_voxel_size['x'])
        prev_voxel_size_y = float(prev_voxel_size['y'])
        prev_voxel_size_z = float(prev_voxel_size['z'])
        nx, ny, nz, nxs, nys, nzs, mx, my, mz = \
            mrc.header.nx, mrc.header.ny, mrc.header.nz, \
            mrc.header.nxstart, mrc.header.nystart, mrc.header.nzstart, \
            mrc.header.mx, mrc.header.my, mrc.header.mz
        orig = mrc.header.origin
        #check the useful density in the input
        input_data = mrc.data
        useful_index = np.argwhere(input_data>contour)
        min_x = int(np.min(useful_index[:,0]))
        max_x = int(np.max(useful_index[:,0]))
        min_y = int(np.min(useful_index[:,1]))
        max_y = int(np.max(useful_index[:,1]))
        min_z = int(np.min(useful_index[:,2]))
        max_z = int(np.max(useful_index[:,2]))
        mapc = mrc.header.mapc
        mapr = mrc.header.mapr
        maps = mrc.header.maps
        new_data = input_data[min_x:max_x,min_y:max_y,min_z:max_z]
        shift_start = permute_map_coord_to_pdb([min_x,min_y,min_z],mapc,mapr,maps)
        origin = np.array(mrc.header.origin.tolist(), dtype=np.float32)
        origin = np.array(origin)+np.array(shift_start)
        mrc_new = mrcfile.new(output_map, data=new_data, overwrite=True)
        vsize = mrc_new.voxel_size
        vsize.flags.writeable = True
        vsize.x = 1.0
        vsize.y = 1.0
        vsize.z = 1.0
        mrc_new.voxel_size = vsize
        mrc_new.update_header_from_data()
        # mrc_new.header.nx = int(max_x-min_x)
        # mrc_new.header.ny = int(max_y-min_y)
        # mrc_new.header.nz = int(max_z-min_z)
        mrc_new.header.nxstart = nxs * prev_voxel_size_x
        mrc_new.header.nystart = nys * prev_voxel_size_y
        mrc_new.header.nzstart = nzs * prev_voxel_size_z
        # mrc_new.header.mx = int(max_x-min_x)
        # mrc_new.header.my = int(max_y-min_y)
        # mrc_new.header.mz = int(max_z-min_z)
        mrc_new.header.mapc = mrc.header.mapc
        mrc_new.header.mapr = mrc.header.mapr
        mrc_new.header.maps = mrc.header.maps
        # mrc_new.header.cella.x = int(max_x-min_x)
        # mrc_new.header.cella.y = int(max_y-min_y)
        # mrc_new.header.cella.z = int(max_z-min_z)
        #mrc_new.header.origin = origin
        (mrc_new.header.origin.x, mrc_new.header.origin.y, mrc_new.header.origin.z) = origin
        mrc_new.update_header_stats()
        mrc.print_header()
        mrc_new.print_header()
        mrc_new.close()
    return output_map
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
