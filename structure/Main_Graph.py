import numpy as np
from structure.Node import Node
from structure.Edge import Edge
import math
import os
from ops.acc_graph import acc_get_density_cut,acc_edge_density,acc_edge_prob
from ops.acc_local_mst import acc_local_mst
class Main_Graph(object):
    def __init__(self,params):
        self.params=params
        self.edge=[]#edge list
        self.Ne=0#number of edges
        self.node=[]#node list
        self.Nnode=0#number of nodes

        self.Nt=0#number of trees
    def copy_gragh(self,graph):
        graph.params = self.params
        graph.edge = self.edge
        graph.Ne = self.Ne
        graph.node = self.node
        graph.Nnode = self.Nnode
        graph.Nt = self.Nt
        graph.cid = self.cid
        graph.adj = self.adj

    def setup_graph(self, point, tmp_data_path):
        """
        :param point: Point class instance
        :param tmp_data_path: path to save intermediate results
        :return:
        """
        self.cid = np.zeros(point.Nmerge)  # chain id list
        self.Nnode = point.Nmerge
        for i in range(self.Nnode):
            tmp_node = Node()
            self.node.append(tmp_node)
        # set the adjunct matrix
        tmp_Nori = int(point.Nori)
        tmp_merged_data = np.array(point.merged_data)
        # tmp_adj=np.array(self.adj)
        tmp_origrid = np.array(point.origrid)
        adj_path = os.path.join(tmp_data_path, 'Quickadj.txt')
        density_path = os.path.join(tmp_data_path, 'Connect_density.txt')
        if os.path.exists(adj_path):
            tmp_adj = np.loadtxt(adj_path)
            tmp_Ne = len(tmp_adj)
        else:
            tmp_adj = np.zeros(self.Nnode * 2000)
            tmp_density = np.zeros(self.Nnode * 2000)
            tmp_Nnode = int(self.Nnode)
            tmp_adj, tmp_Ne, tmp_Density = acc_get_density_cut(tmp_Nori, tmp_merged_data, tmp_adj,
                                                               tmp_density,tmp_origrid, tmp_Nnode)
            Build_connect_dens_dict(tmp_adj, tmp_Density, tmp_data_path)
            tmp_adj = tmp_adj[:tmp_Ne]
            tmp_Density = tmp_Density[:tmp_Ne]
            np.savetxt(density_path, tmp_Density)
            np.savetxt(adj_path, tmp_adj)

        self.adj = set()
        self.Ne = 0
        for i in range(tmp_Ne):
            result = int(tmp_adj[i])
            if result not in self.adj:
                self.adj.add(result)
                i_index = result % self.Nnode
                j_index = int((result - i_index) / self.Nnode)
                if i_index * self.Nnode + j_index not in self.adj:#avoid double adding
                    self.Ne += 1

        print('number of possible connected edges by differen ori data:%d' % self.Ne)

        # Initialize edges
        for i in range(self.Ne):
            tmp_edge = Edge()
            self.edge.append(tmp_edge)
        Ne = 0
        edge_info_path = os.path.join(tmp_data_path, 'alledge_info.txt')
        if not os.path.exists(edge_info_path):
            self.calcu_edge_connect(point, Ne, edge_info_path)
        else:
            for i in range(self.Nnode):
                self.cid[i] = i  # Temparary initialized
            all_edge_info = np.loadtxt(edge_info_path)
            if len(all_edge_info) != self.Ne or (all_edge_info[0, 1] == 0 and all_edge_info[0, 2] == 0):
                all_edge_info = self.calcu_edge_connect(point, Ne, edge_info_path)
            for i in range(self.Ne):
                self.edge[i].d = all_edge_info[i, 0]
                self.edge[i].id1 = int(all_edge_info[i, 1])
                self.edge[i].id2 = int(all_edge_info[i, 2])


    def calcu_edge_connect(self, point, Ne, edge_info_path):
        """
        :param point:
        :param Ne: number of edges
        :param edge_info_path:
        :return:
        """
        for i in range(self.Nnode):
            self.cid[i] = i  # Temparary initialized
            tmp = np.zeros(3)
            for j in range(i + 1, self.Nnode):
                # if self.adj[i,j]==0:
                #    continue
                check_label = i * self.Nnode + j
                if check_label not in self.adj:
                    continue
                d = 0
                for k in range(3):
                    tmp[k] = point.merged_cd_dens[i, k] - point.merged_cd_dens[j, k]
                    d += tmp[k] ** 2
                d = math.sqrt(d)  # two points' distance
                #remove any distance larger than R
                if d>=self.params['R']:
                    print("skip edges with distance %f"%d)
                    continue
                self.edge[Ne].d = d
                self.edge[Ne].id1 = int(i)
                self.edge[Ne].id2 = int(j)
                Ne += 1
        print('sorting edge for MST preparataion %d' % Ne)
        self.Ne = Ne
        #label = self.sort_edge()  # sort edge in ascending order by the distance, prepare for later use
        #if not label:
        #    print('sorting edge did not work, please update code')
        #    return None, None
        all_edge_info = np.zeros([self.Ne, 3])
        for i in range(self.Ne):
            all_edge_info[i, 0] = self.edge[i].d
            all_edge_info[i, 1] = self.edge[i].id1
            all_edge_info[i, 2] = self.edge[i].id2
        np.savetxt(edge_info_path, all_edge_info)
        return all_edge_info

    def sort_edge(self):
        #edge distance to build a list, then sort the list
        distance_list = []
        for i in range(self.Ne):
            tmp_distance = self.edge[i].d
            distance_list.append(tmp_distance)
        distance_list =np.array(distance_list)
        sort_indexes = np.argsort(distance_list)
        final_edge_list = []
        for k in range(len(sort_indexes)):
            cur_index = int(sort_indexes[k])
            final_edge_list.append(self.edge[cur_index])
        self.edge = final_edge_list

        # for i in range(self.Ne):
        #     tmp_edge1 = self.edge[i]
        #     for j in range(i + 1, self.Ne):
        #         tmp_edge2 = self.edge[j]
        #         if tmp_edge1.d > tmp_edge2.d:
        #             tmp_edge = Edge()
        #             tmp_edge.copy_edge(tmp_edge1)
        #             tmp_edge1.copy_edge(tmp_edge2)
        #             tmp_edge2.copy_edge(tmp_edge)
        #             self.edge[j] = tmp_edge2
        #             self.edge[i] = tmp_edge1
        # check the ascending order satisfied or not
        for i in range(self.Ne - 1):
            # if self.edge[i+1].d<self.edge[i].d:
            #    print('%d edge with %f length, %d edge with %f length'%(i,self.edge[i].d,i+1,self.edge[i+1].d))
            #    return False#It makes mistake when 2 d equals same
            if self.edge[i + 1].d == self.edge[i].d:
                self.edge[i + 1].d += 1e-20
        return True

    def set_edge_dens(self,tmp_data_path,mrc,point,density):
        # Set edge density
        # set up filters
        edge_den_path = os.path.join(tmp_data_path, 'edge_dens.txt')
        if not os.path.exists(edge_den_path):
            tmp_dens = self.calcu_edge_dens(mrc, density,point, edge_den_path)
        else:
            tmp_dens = np.loadtxt(edge_den_path)
            edge_den_length = len(tmp_dens)
            if edge_den_length != self.Ne:
                tmp_dens = self.calcu_edge_dens(mrc,density, point, edge_den_path)
        tmp_dens = tmp_dens + 1  # Normalize
        print('edge density min %.2f max %.2f mean %.2f' % (min(tmp_dens), max(tmp_dens), np.mean(tmp_dens)))
        # Here we correct the error of edge dens
        for ii in range(self.Ne):
            self.edge[ii].dens = tmp_dens[ii]
        return tmp_dens

    def calcu_edge_dens(self, mrc,density, point, edge_den_path):
        bandwidth = self.params['g']  # bandwidth of Gaussian distribution
        gstep = mrc.widthx
        fs = (bandwidth / gstep) * 0.5
        fs = fs * fs
        fsiv = 1 / fs
        fmaxd = (bandwidth / gstep) * 2.0
        tmp_Ne = int(self.Ne)
        tmp_id1 = np.zeros(self.Ne)
        tmp_id2 = np.zeros(self.Ne)
        tmp_edge_d = np.zeros(self.Ne)
        for ii in range(self.Ne):
            tmp_id1[ii] = self.edge[ii].id1
            tmp_id2[ii] = self.edge[ii].id2
            tmp_edge_d[ii] = self.edge[ii].d
        tmp_cd_dens = np.array(point.merged_cd_dens)
        ##Notice 09/25/2018,important update, we get normalized density here@@no use any more
        tmp_dens = np.zeros(self.Ne)
        tmp_xdim = int(mrc.xdim)
        tmp_ydim = int(mrc.ydim)
        tmp_zdim = int(mrc.zdim)
        # print(tmp_xdim)
        # exit()
        tmp_mrc_dens = np.array(density)
        tmp_dens = acc_edge_prob(tmp_Ne, tmp_id1, tmp_id2, tmp_cd_dens, tmp_dens, tmp_edge_d, fsiv, fmaxd, tmp_xdim,
                                    tmp_ydim, tmp_zdim, tmp_mrc_dens)
        #tmp_dens = (tmp_dens - np.min(tmp_dens)) / (np.max(tmp_dens) - np.min(tmp_dens)) 

        np.savetxt(edge_den_path, tmp_dens)
        return tmp_dens

    def calcu_local_label(self, point, d2cut):
        cid = np.zeros(point.Nmerge)  # tmp cid for the local MST
        tmp_Nnode = self.Nnode
        tmp_Ne = self.Ne
        tmp_id1 = np.zeros(self.Ne)
        tmp_id2 = np.zeros(self.Ne)
        for ii in range(self.Ne):
            tmp_id1[ii] = self.edge[ii].id1
            tmp_id2[ii] = self.edge[ii].id2
        tmp_cd_dens = np.array(point.merged_cd_dens)
        tmp_local_label = np.zeros(self.Ne)
        tmp_local_label, tmp_count = acc_local_mst(cid, tmp_Nnode, tmp_Ne, tmp_id1, tmp_id2, tmp_cd_dens, d2cut,
                                                   tmp_local_label)
        print('after calulating, with %d local label searches' % tmp_count)
        for ii in range(self.Ne):
            self.edge[ii].local_label = tmp_local_label[ii]

    def config_edge_info(self, point,tmp_data_path):
        dkeep = self.params['k']
        for jj in range(self.Ne):
            self.edge[jj].keep_label = False
            if self.edge[jj].mst_label and self.edge[jj].d <= dkeep:
                self.edge[jj].keep_label = True
        # Init
        for i in range(self.Nnode):
            self.node[i].N = 0
            self.node[i].node_id = i
        tmp_edge_data = np.zeros([self.Ne, 8])  # to save the edge data
        node_path = os.path.join(tmp_data_path, 'node.txt')

        for i in range(self.Ne):
            # Saving data
            tmp_edge_data[i, 0] = self.edge[i].eid
            tmp_edge_data[i, 1] = self.edge[i].id1
            tmp_edge_data[i, 2] = self.edge[i].id2
            tmp_edge_data[i, 3] = self.edge[i].d
            tmp_edge_data[i, 4] = self.edge[i].dens
            tmp_edge_data[i, 5] = int(self.edge[i].mst_label)
            tmp_edge_data[i, 6] = int(self.edge[i].local_label)
            tmp_edge_data[i, 7] = int(self.edge[i].keep_label)
            # initialize the node data
            # if self.edge[i].local_label==False and self.edge[i].mst_label==False:
            #    continue

            id = self.edge[i].id1
            self.node[id].e.append(self.edge[i])

            self.node[id].N += 1
            id = self.edge[i].id2
            self.node[id].e.append(self.edge[i])
            self.node[id].N += 1
            # line=str(self.edge[i].id1)
            # line=line+" "+str(self.edge[i].id2)
            # line=line+" "+str(self.edge[i].eid)
            # print("$$$"+line)
        with open(node_path, 'w') as file:
            for i in range(self.Nnode):
                line = str(self.node[i].node_id) + ' '
                for k in range(4):
                    line = line + str(point.merged_cd_dens[i, k]) + ' '  # x,y,z+density record
                line = line + str(self.node[i].N) + ' '
                for k in range(self.node[i].N):
                    line = line + str(self.node[i].e[k].eid) + ' '
                line = line + '\n'
                file.write(line)
        edge_path = os.path.join(tmp_data_path, 'edge.txt')
        if not os.path.exists(edge_path):
            np.savetxt(edge_path,tmp_edge_data)  # If not, the corruption may happen when some one reading on paralleled situation.
    def assign_node_prob(self,All_Prob_Array,point):
        Name_list = ["atom", "nuc", "chain","base"]
        assert point.Nmerge==self.Nnode
        #use nearby assignment
        ldp_info = point.merged_cd_dens
        check_neighbor = self.params['assign_r']
        for i in range(self.Nnode):
            self.node[i].count_atom = self.params['assign_r']*2+1
            self.node[i].count_nuc = self.params['assign_r']*2+1
            self.node[i].count_chain = self.params['assign_r']*2+1
            self.node[i].count_base = self.params['assign_r']*2+1
            x,y,z,_ = ldp_info[i]
            x, y, z = int(x), int(y), int(z)
            x1, x2 = x-check_neighbor,x+check_neighbor+1
            y1, y2 = y-check_neighbor,y+check_neighbor+1
            z1, z2 = z-check_neighbor,z+check_neighbor+1
            self.node[i].atom_prob = np.sum( All_Prob_Array[Name_list[0]][:,x1:x2,y1:y2,z1:z2],axis=(1,2,3))
            self.node[i].nuc_prob = np.sum( All_Prob_Array[Name_list[1]][:,x1:x2,y1:y2,z1:z2],axis=(1,2,3))
            self.node[i].chain_prob = np.sum(All_Prob_Array[Name_list[2]][:, x1:x2, y1:y2, z1:z2], axis=(1, 2, 3))
            self.node[i].base_prob = np.sum(All_Prob_Array[Name_list[3]][:, x1:x2, y1:y2, z1:z2], axis=(1, 2, 3))
        # init_id, x,y,z,density, merged_to_id
        # merge_info = point.merged_data
        # nouse_count = 0
        # assign_set = set()
        # for i in range(point.Ncd):
        #     merge_to_id,x,y,z,density,merge_status = merge_info[i]
        #     prev_merge_to_id = merge_to_id
        #     use_flag=True
        #     while merge_status==-1:
        #         merge_to_id,_,_,_,_,merge_status = merge_info[int(merge_to_id)]
        #         if merge_to_id == prev_merge_to_id:
        #             use_flag=False
        #             break
        #         prev_merge_to_id = merge_to_id
        #     if use_flag==False or merge_status==-1:
        #         nouse_count +=1
        #         continue
        #     #now the node id: merge_status
        #     x,y,z=int(x),int(y),int(z)
        #     self.node[int(merge_status)].atom_prob += All_Prob_Array[Name_list[0]][:,x,y,z]
        #     self.node[int(merge_status)].count_atom +=1
        #
        #     self.node[int(merge_status)].nuc_prob += All_Prob_Array[Name_list[1]][:, x, y, z]
        #     self.node[int(merge_status)].count_nuc += 1
        #
        #     self.node[int(merge_status)].chain_prob += All_Prob_Array[Name_list[2]][:, x, y, z]
        #     self.node[int(merge_status)].count_chain += 1
        #     assign_set.add(int(merge_status))
        # print("assign prob, we have %d no use ori grid points"%nouse_count)
        # print(assign_set,"we have %d nodes which assigned with prob"%len(assign_set))
        #then check every node to see it's output
        for i in range(self.Nnode):
            # at least one point should be put into the merged node
            assert self.node[i].count_atom!=0
            self.node[i].atom_prob /= self.node[i].count_atom
            self.node[i].nuc_prob /= self.node[i].count_nuc
            self.node[i].chain_prob /= self.node[i].count_chain
            self.node[i].base_prob /= self.node[i].count_base
        sum_patom = np.zeros(5)
        sum_pnuc =  np.zeros(7)
        sum_pchain =  np.zeros(3)
        sum_pbase = np.zeros(5)
        for i in range(self.Nnode):
            sum_patom += self.node[i].atom_prob
            sum_pnuc += self.node[i].nuc_prob
            sum_pchain += self.node[i].chain_prob
            sum_pbase += self.node[i].base_prob
        ref_atom = sum_patom/np.sum(sum_patom)
        ref_nuc = sum_pnuc / np.sum(sum_pnuc)
        ref_chain = sum_pchain / np.sum(sum_pchain)
        ref_base = sum_pbase / np.sum(sum_pbase)
        atom_label_list = ["O", "C", "N", "P", "others"]
        for i in range(len(ref_atom)):
            print("ref %s atom with "%atom_label_list[i],ref_atom[i])
        nuc_label_list = ["sugar", "phosphorus", "A", "UT", "C", "G", "protein"]
        for i in range(len(ref_nuc)):
            print("ref %s nuc with "%nuc_label_list[i],ref_nuc[i])
        chain_label_list = ['DRNA_main', "DRNA_side", "protein"]
        for i in range(len(ref_chain)):
            print("ref %s chain with "%chain_label_list[i],ref_chain[i])
        base_label_list = ["A", "UT", "C", "G", "protein"]
        for i in range(len(ref_base)):
            print("ref %s nuc with " % base_label_list[i], ref_base[i])

        for i in range(self.Nnode):
            self.node[i].log_atom = self.node[i].atom_prob/ref_atom
            self.node[i].log_nuc = self.node[i].nuc_prob / ref_nuc
            self.node[i].log_chain = self.node[i].chain_prob / ref_chain
            self.node[i].log_base = self.node[i].base_prob / ref_base

    def build_subgraph(self,connect_cid):
        unique_id = list(connect_cid)
        unique_id = np.unique(unique_id)
        subgraph =[]
        for select_id in unique_id:
            node_id_list = np.argwhere(connect_cid==select_id)
            if len(node_id_list)<10:
                print("cid %d only included %d nodes,skip!"
                %(select_id,len(node_id_list)))
                continue
            node_id_list = [int(x) for x in node_id_list]
            subgraph.append(list(node_id_list))
        return subgraph
            














import os
def Build_connect_dens_dict(adj_record,density_record,save_path):
    """
    mapping
    :param adj_record:
    :param density_record:
    :param save_path:
    :return:
    a dict mapping [edge_id]:[density]
    """
    save_path=os.path.join(save_path,'Connect_density_dict.txt')
    search_dict={}
    length=len(adj_record)
    for i in range(length):
        search_index=adj_record[i]
        #there are two same way i*Node+j and j*Node+i, therefore,when we use it, we need to check both of them
        if search_index==0:
            continue#no record in adj matrix
        if search_index not in search_dict.keys():
            search_dict[int(search_index)]=density_record[i]
        else:
            search_dict[int(search_index)] = max(density_record[i],search_dict[search_index])
    with open(save_path, 'w') as file:
        file.write(str(search_dict))

    return search_dict


