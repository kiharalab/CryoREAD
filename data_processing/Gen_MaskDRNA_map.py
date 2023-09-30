from data_processing.map_utils import save_dens_map
import mrcfile
import numpy as np
def Gen_MaskDRNA_map(chain_prob,cur_map_path,save_map_path,contour,threshold=0.6):
    sp_prob = chain_prob[0]+chain_prob[1]
    base_prob = chain_prob[-1]
    protein_prob = chain_prob[-2]
    with mrcfile.open(cur_map_path,permissive=True) as mrc:
        dens_data=np.array(mrc.data)
    dens_data[protein_prob>=threshold]=0
    drna_predictions= sp_prob+base_prob
    dens_data[drna_predictions<=0.1]=0
    dens_data[dens_data<=contour]=0
    #then save the new density data
    save_dens_map(save_map_path,dens_data, cur_map_path)


def Gen_MaskProtein_map(chain_prob,cur_map_path,save_map_path,contour,threshold=0.3):
    protein_prob = chain_prob[-2]
    with mrcfile.open(cur_map_path,permissive=True) as mrc:
        dens_data=np.array(mrc.data)
    dens_data[protein_prob<=threshold]=0
    dens_data[dens_data<=contour]=0
    #then save the new density data
    save_dens_map(save_map_path,dens_data, cur_map_path)
