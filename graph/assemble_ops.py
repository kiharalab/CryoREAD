

import numpy as np

import os

def bool_intersection(interval1,interval2):
    ldp_intersection_starting_index = max(interval1[0],interval2[0])
    ldp_intersection_ending_index = min(interval1[1],interval2[1])
    if ldp_intersection_ending_index>ldp_intersection_starting_index:
        return True
    else:
        return False

def overlap_size(interval1,interval2):
    ldp_intersection_starting_index = max(interval1[0],interval2[0])
    ldp_intersection_ending_index = min(interval1[1],interval2[1])
    if ldp_intersection_ending_index>ldp_intersection_starting_index:
        return ldp_intersection_ending_index-ldp_intersection_starting_index
    else:
        return 0


def calculate_identity(match_align1_seq,match_align2_seq):
    count=0
    assert len(match_align1_seq)==len(match_align2_seq)
    for k in range(len(match_align1_seq)):
        if match_align1_seq[k]==match_align2_seq[k]:
            count+=1
    return count/len(match_align1_seq)
def build_collision_table(All_Base_Assign_List,checking_stride,ldp_size,
            overall_dict,collision_save_path,soft_rule=1,identity_cutoff=0.8,
            load_collision=True,overlap_cutoff=0.2):
    #build n*n collision table for the assembling
    order_key_index={}#map current order to original key
    order_chain_index= {}#to index which chain it corressponds to
    key_order_index ={} #map original key to current order
    order_index = 0
    count_total_frag=0
    #1. collect all relationship of chain
    to_location_idx ={}
    for k,cur_path_list in enumerate(All_Base_Assign_List):
        current_base_list = All_Base_Assign_List[k]
        start_ldp_index =0
        while start_ldp_index<len(current_base_list):
            current_key = "%d_%d"%(k,start_ldp_index)
            final_path_assign_collection_list = overall_dict[current_key]

            for j in range(len(final_path_assign_collection_list)):
                order_key_index[order_index]=current_key
                key_order_index[current_key]= order_index
                order_chain_index[order_index]=j
                to_location_idx[order_index]=count_total_frag
                order_index+=1

            start_ldp_index+=checking_stride

            count_total_frag+=1

    N_fragement = order_index
    if os.path.exists(collision_save_path):
        if load_collision:
            collision_table=np.load(collision_save_path)
        else:
            collision_table=None
        if collision_table is None:
            return None,order_key_index,order_chain_index,key_order_index
        if len(collision_table)==N_fragement:
            return collision_table,order_key_index,order_chain_index,key_order_index
    #2. build collision table based on overlap checking.
    print("in total we have %d fragements to assemble"%N_fragement)
    collision_table = np.zeros([N_fragement,N_fragement]).astype(int)#1 denotes forbidding co-existance
    count_ldp_collision1 = 0
    count_ldp_collision2 = 0
    count_ldp_collision3 = 0
    count_ldp_collision4 = 0
    count_ldp_acceptable =0
    count_seq_collision1=0
    count_seq_collision2=0
    count_seq_collision3=0
    count_seq_collision4=0
    count_seq_acceptable=0
    for i in range(N_fragement):
        current_key1 = order_key_index[i]
        current_chain_candidate1 = order_chain_index[i]
        split_key1 = current_key1.split("_")
        ldp_path1= int(split_key1[0])
        ldp_starting_index1 = int(split_key1[1])
        ldp_ending_index1 = ldp_starting_index1+ldp_size
        current_seq_info1 = overall_dict[current_key1][current_chain_candidate1]
        location_idx1 = to_location_idx[i]

        for j in range(1,N_fragement):
            current_key2 = order_key_index[j]
            current_chain_candidate2 = order_chain_index[j]
            split_key2 = current_key2.split("_")
            ldp_path2= int(split_key2[0])
            ldp_starting_index2 = int(split_key2[1])
            ldp_ending_index2 = ldp_starting_index2+ldp_size
            current_seq_info2 = overall_dict[current_key2][current_chain_candidate2]
            location_idx2 = to_location_idx[j]
            ldp_overlap=False
            #check if any overlap between two ldps
            if ldp_path1==ldp_path2:#if they are in the same path tracing
                ldp_union_starting_index = min(ldp_starting_index1,ldp_starting_index2)
                ldp_union_ending_index = max(ldp_ending_index1,ldp_ending_index2)
                if ldp_union_ending_index-ldp_union_starting_index<2*ldp_size:
                    #path has overlap, we need to mark collisions
                    ldp_overlap=True
                    ldp_intersection_start = max(ldp_starting_index1,ldp_starting_index2)
                    ldp_intersection_end =  min(ldp_ending_index1,ldp_ending_index2)
                    if ldp_union_starting_index==ldp_starting_index1:
                        match_info1 = current_seq_info1
                        match_info2 = current_seq_info2
                    else:
                        match_info1 = current_seq_info2
                        match_info2 = current_seq_info1
                    match_intersect_seq1 = match_info1['match_seq']
                    match_intersect_seq2 = match_info2['match_seq']
                    match_intersect_chain1 = match_info1['chain']
                    match_intersect_chain2 = match_info2['chain']
                    match_intersect_direction1 = match_info1['direction']
                    match_intersect_direction2 = match_info2['direction']
                    match_intersect_interval1 = match_info1['interval']
                    match_intersect_interval2 = match_info2['interval']
                    #2.1 overlaped ldp regions should have same assignment direction
                    if match_intersect_direction1!=match_intersect_direction2 or match_intersect_chain1!=match_intersect_chain2:
                        collision_table[i,j]=1
                        collision_table[j,i]=1
                        count_ldp_collision1+=1
                        continue
                    #2.2 check the if there exist assignment intersection or not, they must have intersections
                    if not bool_intersection(match_intersect_interval1,match_intersect_interval2):
                        collision_table[i,j]=1
                        collision_table[j,i]=1
                        count_ldp_collision2+=1
                        continue
                    #if same chain, then check the interval overlap to compare the intersection
                    match_align1_begin = ldp_intersection_start-ldp_union_starting_index
                    match_align1_useful_count = len([k for k in match_intersect_seq1[:match_align1_begin] if k!="-"])
                    #calculate the sequence assignment of match 1 aligned interval
                    match_align1_interval = [match_intersect_interval1[0]+match_align1_useful_count,match_intersect_interval1[1]]
                    match_align1_seq = match_intersect_seq1[match_align1_begin:]

                    match_align2_end = ldp_intersection_end-ldp_intersection_start
                    match_align2_useful_count = len([k for k in match_intersect_seq2[:match_align2_end] if k!="-"])
                    match_align2_interval = [match_intersect_interval2[0],match_intersect_interval2[0]+match_align2_useful_count]
                    match_align2_seq = match_intersect_seq2[:match_align2_end]
                    #2.3 for overlap regions, must be overlap in sequence also check the if intersection or not, they must have intersections
                    if not bool_intersection(match_align1_interval,match_align2_interval):
                        collision_table[i,j]=1
                        collision_table[j,i]=1
                        count_ldp_collision3+=1
                        continue
                    #check the overlap sequence to see if they can agree with each other
                    cur_identity=calculate_identity(match_align1_seq,match_align2_seq)
                    if cur_identity<identity_cutoff:
                        collision_table[i,j]=1
                        collision_table[j,i]=1
                        count_ldp_collision4+=1
                        continue
                    count_ldp_acceptable+=1
            if ldp_overlap is False and soft_rule==0:
                #no overlap in ldp fragment, then they must not from the same sequence
                match_info1 = current_seq_info1
                match_info2 = current_seq_info2
                match_intersect_chain1 = match_info1['chain']
                match_intersect_chain2 = match_info2['chain']
                match_intersect_direction1 = match_info1['direction']
                match_intersect_direction2 = match_info2['direction']
                match_intersect_interval1 = match_info1['interval']
                match_intersect_interval2 = match_info2['interval']
                match_seqfrag1_length = abs(match_intersect_interval1[1]-match_intersect_interval1[0])
                match_seqfrag2_length = abs(match_intersect_interval2[1]-match_intersect_interval2[0])
                max_seqfrag_length = max(match_seqfrag1_length,match_seqfrag2_length)
                match_intersect_chain1length = match_info1['chain_length']
                match_intersect_chain2length = match_info2['chain_length']
                if match_intersect_chain1!=match_intersect_chain2:
                    #if they are from different chain, totally fine
                    count_seq_acceptable+=1
                    continue
                #now it's the same chain
                #if it's reverse direction, we need to check if they have overlap or not
                #if overlap,we only favor 1
                if match_intersect_direction1!=match_intersect_direction2:
                    match_chain1_new_interval = [match_intersect_chain1length-match_intersect_interval1[1]-1,match_intersect_chain1length-match_intersect_interval1[0]]
                    if overlap_size(match_chain1_new_interval, match_intersect_interval2)>=max_seqfrag_length*overlap_cutoff:
                        collision_table[i,j]=1
                        collision_table[j,i]=1
                        count_seq_collision1+=1
                        continue
                else:
                    #if it's same direction and have overlap, then avoid assign
                    if overlap_size(match_intersect_interval1, match_intersect_interval2)>=max_seqfrag_length*overlap_cutoff:
                        collision_table[i,j]=1
                        collision_table[j,i]=1
                        count_seq_collision2 +=1

                        continue

                count_seq_acceptable+=1
        if i%10==0:
            print("collision table finished %d/%d"%(i,len(collision_table)))
    print("-"*30)
    print("ldp overlap collision report:")
    print("1. different chain or different direction:%d"%count_ldp_collision1)
    print("2. no overlap in sequence:%d"%count_ldp_collision2)
    print("3. matched region sequence no overlap:%d"%count_ldp_collision3)
    print("4. matched region sequence low identity:%d"%count_ldp_collision4)
    print("# Total acceptable pass ldp overlap check: %d"%count_ldp_acceptable)

    print("-"*30)
    print("seq overlap collision report:")
    print("1. interaction of different direction assignment:%d"%count_seq_collision1)
    print("2. interaction of same direction assignment:%d"%count_seq_collision2)
    print("3. nearby assignment did not satisfy distance constraint:%d"%count_seq_collision3)
    print("# Total acceptable pass for seq in same chain:%d"%count_seq_acceptable)
    print("-"*30)
    np.save(collision_save_path,collision_table)
    return collision_table,order_key_index,order_chain_index,key_order_index


from ortools.sat.python import cp_model
from ortools.linear_solver import pywraplp
#from ortools.init import pywrapinit
import numpy as np
def prepare_score_array(N_order, order_key_index,order_chain_index,overall_dict):
    score_array = np.zeros(N_order)
    for k in range(N_order):
        current_key1 = order_key_index[k]
        current_chain_candidate1 = order_chain_index[k]
        seq_info = overall_dict[current_key1][current_chain_candidate1]
        current_score = seq_info['score']
        score_array[k]=int(current_score*100)#keep 0.01 precision
    return score_array

import time 
def solve_assignment(collision_table,order_key_index,order_chain_index,overall_dict,time_use=3600):
    model = cp_model.CpModel()

    x = []
    for i in range(len(collision_table)):
        x.append(model.NewBoolVar(f'x[{i}]'))

    print('#Adding constraints')
    #collision can't  be true, if current assignment is determined
    for i in range(len(collision_table)):
        model.Add(sum(x[j] for j in range(len(collision_table)) if collision_table[i][j]==1)==0).OnlyEnforceIf(x[i])

    score_array=prepare_score_array(len(collision_table), order_key_index,order_chain_index,overall_dict)
    score_array = score_array.astype(int)
    print('#Setting Objective Score')
    # Objective
    objective_terms = []
    #[fid,rawsco,zsco,ali,tabu]
    for i in range(len(collision_table)):
        objective_terms.append(score_array[i]* x[i]) #Sum of Raw Scores
    model.Maximize(sum(objective_terms))


    print('#Start Solving...')
    # Solve
    print("Time limit is set to %d seconds"%time_use)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_use
    solver.parameters.num_search_workers = 8
    solver.parameters.log_search_progress = True
    solution_printer = cp_model.ObjectiveSolutionPrinter()
    time1 = time.time()
    status = solver.SolveWithSolutionCallback(model, solution_printer)
    time2 = time.time()
    print("Time used for solving:",time2-time1)
    """
    %unignore operations_research::MPSolver::ResultStatus;
    %unignore operations_research::MPSolver::OPTIMAL; value 0
    %unignore operations_research::MPSolver::FEASIBLE;  value 1 // No unit test
    %unignore operations_research::MPSolver::INFEASIBLE; value 2
    %unignore operations_research::MPSolver::UNBOUNDED; value 3 // No unit test
    %unignore operations_research::MPSolver::ABNORMAL; value 4
    %unignore operations_research::MPSolver::NOT_SOLVED; value 5 // No unit test

    OPTIMAL = _pywraplp.Solver_OPTIMAL
    r optimal
    FEASIBLE = _pywraplp.Solver_FEASIBLE
    r feasible, or stopped by limit.
    INFEASIBLE = _pywraplp.Solver_INFEASIBLE
    r proven infeasible.
    UNBOUNDED = _pywraplp.Solver_UNBOUNDED
    r proven unbounded.
    ABNORMAL = _pywraplp.Solver_ABNORMAL
    r abnormal, i.e., error of some kind.
    NOT_SOLVED = _pywraplp.Solver_NOT_SOLVED
    r not been solved yet.
    """
    print("current status:",status)
    results=[]
    if status == pywraplp.Solver.OPTIMAL or status==pywraplp.Solver.FEASIBLE:
        try:
            for i in range(len(collision_table)):
                if solver.BooleanValue(x[i]):
                    #print(i,'ID=',idtbl[i]) #order->ID
                    results.append(i)
        except:
            print("solution status is optimal or feasible, but error raised. ",results)
            results=[]
            
    # elif status == pywraplp.Solver.INFEASIBLE  or status == pywraplp.Solver.UNBOUNDED:
    #     print("current status:",status)
    #     print("no optimal or feasible solution found for assembling, use temporary solution for final results.")
    #     try:
    #         for i in range(len(collision_table)):
    #             if solver.BooleanValue(x[i]):
    #                 #print(i,'ID=',idtbl[i]) #order->ID
    #                 results.append(i)
    #     except:
    #         return []#indicate falure
    else:
        print("*"*100)
        print("current status:",status)
        print('No solution found for assembling, please make contact the developer! You can also check the CryoREAD_noseq.pdb as temporary results.')
        print("*"*100)
    return results
