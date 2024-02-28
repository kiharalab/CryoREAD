#make this file self consistent without too many dependency
from collections import defaultdict
import numpy as np
from scipy.spatial.distance import cdist
def cif2dict(input_cif_path,filter_list=None):      
    """
    input_cif_path: input cif file path
    return:
    a dictionary in this format: [nuc_id][atom_id]:[coordinates]
    """
    begin_check=False
    block_list=[]
    with open(input_cif_path,'r') as rfile:
        for line in rfile:
            if "loop_" in line:
                begin_check=True
                continue

            if begin_check and "_atom_site" in line:
                block_list.append(line.strip("\n").replace(" ",""))
                continue
            if begin_check and "_atom_site" not in line:
                begin_check=False
    atom_ids = block_list.index('_atom_site.id')
    try:
        seq_ids = block_list.index('_atom_site.label_seq_id')
    except:
        seq_ids = block_list.index('_atom_site.auth_seq_id')
    try:
        chain_ids = block_list.index("_atom_site.auth_asym_id")
    except:
        chain_ids = block_list.index("_atom_site.label_asym_id")
    atom_type_ids = block_list.index("_atom_site.label_atom_id")
    res_name_ids = block_list.index("_atom_site.label_comp_id")
    x_ids = block_list.index("_atom_site.Cartn_x")
    y_ids = block_list.index("_atom_site.Cartn_y")
    z_ids = block_list.index("_atom_site.Cartn_z")
    structure_dict =defaultdict(dict)#[chain_id][nuc_id]:[coordinates, nuc_type,nuc_score]
   
    map_dict={"A":0,"U":1,"T":1,"C":2,"G":3,"DA":0,"DU":1,"DT":1,"DC":2,"DG":3}
    total_nuc_id=0
    track_nuc_id=-1000
    with open(input_cif_path,'r') as rfile:
        for line in rfile:
            if line.startswith("ATOM"):
                split_result = line.split()
                split_info=line.strip("\n").split()
                current_atom_name = split_info[atom_type_ids]
                current_atom_name = current_atom_name.replace(" ","")
                current_res_index = int(split_info[seq_ids])
                current_res_name = split_info[res_name_ids]
                current_x = float(split_info[x_ids])
                current_y = float(split_info[y_ids])
                current_z = float(split_info[z_ids])
                current_atom_name = current_atom_name.replace(" ","")
                if current_res_name not in map_dict:
                    continue
                pred_label = map_dict[current_res_name]
                
                if filter_list is not None and current_atom_name not in filter_list:
                    continue
                if track_nuc_id!=current_res_index:
                    track_nuc_id=current_res_index
                    total_nuc_id+=1
                
                structure_dict[total_nuc_id][current_atom_name]=[current_x,current_y,current_z,pred_label]
                
                
        return structure_dict

def pdb2dict(input_pdb_path,filter_list=None):
    """
    input_pdb_path: input pdb file path
    return:
    a dictionary in this format: [nuc_id][atom_id]:[coordinates]
    
    """
    structure_dict =defaultdict(dict)#[chain_id][nuc_id]:[coordinates, nuc_type,nuc_score]
    map_dict={"A":0,"U":1,"T":1,"C":2,"G":3,"DA":0,"DU":1,"DT":1,"DC":2,"DG":3}
    total_nuc_id=0
    track_nuc_id=-1000
    with open(input_pdb_path,'r') as rfile:
        for line in rfile:
            if line.startswith("ATOM"):
                chain_id = line[21]
                atom_name = line[12:16]
                x=float(line[30:38])
                y=float(line[38:46])
                z=float(line[46:54])
                nuc_id=int(line[22:26])
                resn = line[17:20]
                coordinates = [x,y,z]
                nuc_type = resn#split_result[5]
                nuc_type = nuc_type.replace(" ","")
                if nuc_type not in map_dict:
                    continue
                pred_label = map_dict[nuc_type]
                atom_name = atom_name.replace(" ","")
                if filter_list is not None and atom_name not in filter_list:
                    continue
                if track_nuc_id!=nuc_id:
                    track_nuc_id=nuc_id
                    total_nuc_id+=1
                
                structure_dict[total_nuc_id][atom_name]=[x,y,z,pred_label]
                
                
        return structure_dict
def calcudate_atomwise_distmat(query_dict,target_dict):
    """
    query_dict: dictionary of predicted structure
    target_dict: dictionary of native structure
    return:
    distance_matrix: M*N matrix, M is the number of nucleotides in query pdb, N is the number of nucleotides in target pdb
    """
    #first get atom list
    query_keys=list(query_dict.keys())
    target_keys = list(target_dict.keys())
    tmp_key = query_keys[0]
    atom_list = list(target_dict[tmp_key].keys())
    distance_matrix = np.zeros([len(query_keys),len(target_keys)]) #M*N matrix
    count_matrix = np.zeros([len(query_keys),len(target_keys)]) #M*N matrix count the number of atoms for each pair of nucleotides
    for select_atom in atom_list:
        query_atom_list=[]
        missing_query_index=[]
        for k,nuc_id in enumerate(query_keys):
            if select_atom not in query_dict[nuc_id]:
                
                print("atom %s not found in target pdb"%select_atom)
                query_atom_list.append([999999,999999,999999])#indicate record missing atom
                missing_query_index.append(k)
            else:
                query_atom_list.append(query_dict[nuc_id][select_atom][:3])
        target_atom_list=[]
        missing_target_index=[]
        for k,nuc_id in enumerate(target_keys):
            if select_atom not in target_dict[nuc_id]:
                print("atom %s not found in query pdb"%select_atom)
                target_atom_list.append([999999,999999,999999])
                missing_target_index.append(k)
            else:
                target_atom_list.append(target_dict[nuc_id][select_atom][:3])
        query_atom_list=np.array(query_atom_list)
        target_atom_list=np.array(target_atom_list)
        current_dist_matrix = cdist(query_atom_list,target_atom_list) #M*N distance matrix for current atom type
        for k in range(len(missing_query_index)):
            current_dist_matrix[missing_query_index[k],:]=0
        for k in range(len(missing_target_index)):
            current_dist_matrix[:,missing_target_index[k]]=0
        distance_matrix+=current_dist_matrix
        tmp_count_matrix = np.ones([len(query_keys),len(target_keys)])
        for k in range(len(missing_query_index)):
            tmp_count_matrix[missing_query_index[k],:]=0
        for k in range(len(missing_target_index)):  
            tmp_count_matrix[:,missing_target_index[k]]=0
        count_matrix+=tmp_count_matrix
    #check zero in count matrix
    zero_index = np.where(count_matrix==0)
    if len(zero_index[0])>0:
        print("Some atoms are missing in the query/target pdb/cif, please check the input file")
        exit()
    distance_matrix = distance_matrix/count_matrix
    return distance_matrix

def calculate_eval_metric(query_dict,target_dict,distance_matrix,cutoff):
    """
    query_dict: dictionary of predicted structure
    target_dict: dictionary of native structure 
    distance_matrix: M*N matrix, M is the number of nucleotides in query pdb, N is the number of nucleotides in target pdb
    cutoff: distance cutoff for evaluation
    return:
    atom_coverage,atom_precision,sequence_match,sequence_match_prec,sequence_recall,sequence_prec,RMSD
    """
    

    #first find the closest pair of nucleotides
    query_keys=list(query_dict.keys())
    target_keys = list(target_dict.keys())
    #calculate atom coverage if we simply find the closest pair of nucleotides
    visit_target_set=set()
    visit_query_set=set()
    all_possible_match_dict = defaultdict(list)
    RMSD=[]
    for j,target_key in enumerate(target_keys):
        cur_dist = distance_matrix[:,j]
        current_query_index = np.argwhere(cur_dist<=cutoff)

        if len(current_query_index)>0:
            for tmp_index in current_query_index:
                tmp_index = int(tmp_index)
                visit_query_set.add(tmp_index)
                all_possible_match_dict[j].append(tmp_index)
            visit_target_set.add(j)
            current_query_index = int(np.argmin(cur_dist))
            RMSD.append(np.min(cur_dist))
    RMSD= np.mean(np.array(RMSD))
    visit_atom_query_set = visit_query_set
    visit_atom_target_set = visit_target_set
    atom_precision = len(visit_query_set)/len(query_keys)
    atom_coverage = len(visit_target_set)/len(target_keys)
    #get the atom coverage and precision
    overall_match_dict={} #match_dict[query_nuc_id]=target_nuc_id
    visit_query_set=set()#indicate if the target is already matched
    visit_target_set=set()
    #first check the closest pair of nucleotides
    for j,target_key in enumerate(target_keys):
        cur_dist = distance_matrix[:,j]
        current_query_index = int(np.argmin(cur_dist))
        query_key = query_keys[current_query_index]
        match_dist = cur_dist[current_query_index]
        if match_dist>cutoff:
            continue   
        if current_query_index in visit_query_set or j in visit_target_set:
            continue
        current_query_base = query_dict[query_key]["P"][3]
        current_target_base = target_dict[target_key]["P"][3]
        if current_query_base!=current_target_base:
            continue
        #1 predict base can only be mapped to one native base
        overall_match_dict[query_key]=target_key
        visit_query_set.add(current_query_index)
        visit_target_set.add(j)
    #then check all the remained nucleotides in target that did not find match
    for j,target_key in enumerate(target_keys):
        if j in visit_target_set:
            continue
        cur_dist = distance_matrix[:,j]
       
        current_target_base = target_dict[target_key]["P"][3]
        current_query_index = np.argwhere(cur_dist<=cutoff)
        if len(current_query_index)>0:
            for tmp_index in current_query_index:
                tmp_index = int(tmp_index)
                if tmp_index in visit_query_set:
                    continue
                query_key = query_keys[tmp_index]
                current_query_base = query_dict[query_key]["P"][3]
                if current_query_base==current_target_base:
                    overall_match_dict[query_key]=target_key
                    visit_query_set.add(tmp_index)
                    visit_target_set.add(j)
                    break
    total_match = len(overall_match_dict)
    sequence_match = total_match/len(visit_atom_target_set)
    sequence_match_prec = total_match/len(visit_atom_query_set)
    sequence_recall = total_match/len(target_keys)
    sequence_prec = total_match/len(query_keys)
    #calculate rmsd based on the matched nucleotides
    return atom_coverage,atom_precision,sequence_match,sequence_match_prec,sequence_recall,sequence_prec,RMSD

def evaluate_structure(query_pdb,target_pdb,cutoff=5.0):
    
    """
    query_pdb: predicted pdb file
    target_pdb: native pdb file
    cutoff: distance cutoff for evaluation
    """
    _pho_atoms=["OP3", "P", "OP1", "OP2", "O5'",]
    _sugar_atoms = [ "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'"]    

    filter_list = _pho_atoms+_sugar_atoms+["N1"] #backbone atoms for evaluation
    if query_pdb.endswith(".cif"):
        query_dict = cif2dict(query_pdb,filter_list)
    elif query_pdb.endswith(".pdb"):
        query_dict = pdb2dict(query_pdb,filter_list)
    else:
        print("format of query file %s not supported"%query_pdb)
        exit()
    print("Input query pdb includes %d nucleotides"%len(query_dict))
    if len(query_dict)==0:
        print("No nucleotide found in the query pdb/cif file")
        exit()
    #read the target pdb file
    if target_pdb.endswith(".cif"):
        target_dict = cif2dict(target_pdb,filter_list)
    elif target_pdb.endswith(".pdb"):
        target_dict = pdb2dict(target_pdb,filter_list)
    else:
        print("format of target file %s not supported"%target_pdb)
        exit()
    print("Input target pdb includes %d nucleotides"%len(target_dict))
    if len(target_dict)==0:
        print("No nucleotide found in the target pdb/cif file")
        exit()
    #calculate the distance between the two structures: N*M distance matrix, this is the average distance of corresponding atoms between the two structures 
    distance_matrix = calcudate_atomwise_distmat(query_dict,target_dict)
    #calculate atom coverage, atom precision, sequuence recall(match),sequence precision(match), sequence recall, sequence precision, RMSD, base-RMSD
    atom_coverage,atom_precision,sequence_match,sequence_match_prec,sequence_recall,sequence_prec,RMSD = calculate_eval_metric(query_dict,target_dict,distance_matrix,cutoff)
    print("*"*100)
    print("Atom Coverage: %.3f"%atom_coverage+" Atom Precision: %.3f"%atom_precision)
    print("Sequence Recall(Match): %.3f"%sequence_match+" Sequence Precision(Match): %.3f"%sequence_match_prec)
    print("Sequence Recall: %.3f"%sequence_recall+" Sequence Precision: %.3f"%sequence_prec)
    print("RMSD: %.3f"%RMSD)
    print("*"*100)