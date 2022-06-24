
from data_processing.map_utils import process_map_data
import numpy as np
from graph.MRC import MRC
from graph.LDP_ops import build_LDP

def Build_Unet_Graph(origin_map_path,chain_prob_path,fasta_path,save_path,
                     gaussian_bandwidth,dcut,rdcut, params):
    pho_prob_threshold=0.1#make sure all necessary edges are connected
    sugar_prob_threshold = 0.1
    base_prob_threshold = 0.25
    map_data, mapc, mapr, maps, origin, nxstart, nystart, nzstart = process_map_data(origin_map_path)
    map_info_list=[mapc, mapr, maps, origin, nxstart, nystart, nzstart]
    chain_class =8

    chain_prob = np.load(chain_prob_path)#[sugar,phosphate,A,UT,C,G,protein,base]
    input_mrc = MRC(origin_map_path, gaussian_bandwidth)

    #1. chain tracing
    sp_prob = chain_prob[0]+chain_prob[1]
    pho_prob = chain_prob[1]
    sugar_prob = chain_prob[0]

    input_mrc.upsampling_pho_prob(pho_prob,threshold=pho_prob_threshold,filter_array=None)
    input_mrc.upsampling_sugar_prob(sugar_prob,threshold=sugar_prob_threshold,filter_array=None)

    #1.1 LDP construction based on probability map
    pho_point= build_LDP(input_mrc,input_mrc.pho_dens, input_mrc.pho_Nact,origin_map_path,save_path,"pho",params,map_info_list)
    sugar_point = build_LDP(input_mrc,input_mrc.sugar_dens,input_mrc.sugar_Nact,origin_map_path,save_path,"sugar",params,map_info_list)

    
