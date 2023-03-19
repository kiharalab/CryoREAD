
import numpy as np
from structure.Node import Node
import os
class Tree(object):
    def __init__(self,params):
        self.params=params
        self.len=0#tree length
        self.bf_len=0
        self.Nnode=0
        self.node=[]#node list
        self.Ntotal=0#number of nodes
        self.Etotal=0#number of edges
        self.Ne=0#Number of edges in tree
        self.St=0#Start point
        self.Ed=0
        #self.ActE=[]#Active edge list:label True or False
        #self.ActN=[]#Active node list:label True or False
        self.Nstock=0
        self.Lpath=0#path length
        self.Nadd=0
        self.Ncut=0
        self.score=0
        self.Nmv=0#Records the move list length

    def Inittree(self,label):
        if label:
            if self.Nnode==0 or self.Etotal==0 or self.Ne==0:
                print("Not necessary init for the tree")
                return True


            self.AddTbl=np.zeros(self.Etotal)
            self.CutTbl=np.zeros(self.Etotal)
            self.stock=np.zeros(self.Nnode)

            self.cid=np.zeros(self.Nnode)
            self.cost=np.zeros(self.Nnode)
            #self.node=[]
            self.mv=[]
            #for i in range(self.Ne**2):
            #    tmp_Move=Move()
            #    self.mv.append(tmp_Move)
            self.Path=np.zeros(self.Ne)

            for i in range(self.Nnode):
                tmp_node=Node()
                tmp_node.N=0
                self.node.append(tmp_node)
        self.ActN=np.zeros(self.Nnode)
        self.nextv=np.zeros(self.Nnode)
        self.cost=np.zeros(self.Nnode)
        self.ActE=np.zeros(self.Etotal)
        self.MaskN = np.zeros(self.Nnode)
        self.UsedN = np.zeros(self.Nnode)
        self.Path = np.zeros(self.Ne+1)
        return False
    def CopyTree(self,tree_in,label):
        self.len=tree_in.len
        self.Nnode=tree_in.Nnode
        self.Ne=tree_in.Ne
        self.Ntotal=tree_in.Ntotal
        self.Etotal=tree_in.Etotal
        self.St=tree_in.St
        self.Ed=tree_in.Ed
        self.Nadd=tree_in.Nadd
        self.Ncut=tree_in.Ncut
        self.score=tree_in.score
        self.Lpath=tree_in.Lpath
        self.Nmv=tree_in.Nmv
        self.ActE=np.zeros(self.Etotal)
        self.MaskN = np.zeros(tree_in.Nnode)
        self.UsedN = np.zeros(tree_in.Nnode)
        if label==True:
            #self.node=[]
            #for i in range(self.Nnode):
            #    tmp_node=Node()
            #    self.node.append(tmp_node)

            self.AddTbl=np.zeros(self.Etotal)
            self.CutTbl=np.zeros(self.Etotal)
            self.stock=np.zeros(self.Nnode)

            self.cid=np.zeros(self.Nnode)

            self.node=[]
            self.mv=[]
            #for i in range(self.Etotal**2):
            #    tmp_Move=Move()
            #    self.mv.append(tmp_Move)
        self.cost=np.zeros(self.Nnode)
        self.ActN=np.zeros(self.Nnode)
        self.Path=np.zeros(self.Ne)
        self.nextv=np.zeros(self.Nnode)
        for i in range(self.Lpath):
            self.Path[i]=tree_in.Path[i]
        for i in range(self.Etotal):
            self.ActE[i]=tree_in.ActE[i]
        for i in range(self.Nnode):
            self.cost[i]=tree_in.cost[i]
            self.ActN[i]=tree_in.ActN[i]
            self.nextv[i]=tree_in.nextv[i]
    def MoveTree(self,graph,cut_id,add_id):

        #Update len
        self.len=self.len - graph.edge[cut_id].d + graph.edge[add_id].d
        #Update ActiveE
        self.ActE[cut_id]=0
        self.ActE[add_id]=1

    def Setup_Connection(self,graph):
        Nmin_ldp = 10
        label = graph.sort_edge()  # sort edge in ascending order by the distance
        if not label:
            print('1st sorting edge did not work, please update code')
            return True
            # Clean trees

        MaxCid = 0
        Nt = 0  # id for tree
        tree = []
        print('finding mst')
        for i in range(graph.Ne):
            v1 = graph.edge[i].id1
            v2 = graph.edge[i].id2
            if graph.cid[v1] == graph.cid[v2]:
                #they already put into the connected region
                continue
            tree.append(graph.edge[i])
            tmp_cid = graph.cid[v2]
            Nt += 1
            if MaxCid < tmp_cid:
                MaxCid = tmp_cid
            # if smaller edge shares one point with another larger edge,
            # if one point of the smaller edge's chain id equals to the larger id point of the larger edge,
            #  we change this smaller edge's point's chain id to the smaller id point of the large edge
            # You can understand it as the connection id,which smaller distance 1 connection id connects to the larger one another point
            # Chain id is here to avoid loop in mininmum tree
            for j in range(Nt):
                #current edge is also in the checking loop
                if graph.cid[tree[j].id1] == tmp_cid:
                    graph.cid[tree[j].id1] = graph.cid[v1]  # It can only use we keep id1<id2
                if graph.cid[tree[j].id2] == tmp_cid:
                    graph.cid[tree[j].id2] = graph.cid[v1]
        graph.tree = tree
        graph.Nt = Nt
        print('Found the mst wanted')
        print('Maxcid %d' % MaxCid)
        print('Found now the Nt is %d' % (graph.Nt))
        Ncid = np.zeros(int(MaxCid) + 1)
        for i in range(graph.Nnode):
            if graph.cid[i] <= MaxCid:
                Ncid[int(graph.cid[i])] += 1  # Counts the points how many connections are, including the connection through 2 or 3 edges
        UseCid = -1  # Use cid is to record the self.cid[i] with most connections
        Nuse = 0
        for i in range(int(MaxCid) + 1):
            if Nuse < Ncid[i]:
                UseCid = i
                Nuse = Ncid[i]  # Nuse is to find the largest
            if Ncid[i]>0:
                print("&#CID %d N= %d\n"%(i,Ncid[i]))
        if UseCid == -1: # this is the cluster that will be used to construct mst
            return True  # Break the running process

        Ntmp = 0
        for i in range(graph.Ne):
            if Ncid[int(graph.cid[graph.edge[i].id1])]>=Nmin_ldp:
                graph.edge[Ntmp] = graph.edge[i]
                Ntmp += 1
        graph.Ne = Ntmp  # tree edges
        # Sort the edge
        graph.edge = graph.edge[0:graph.Ne]  # Only keep the edges in the graph

        print('finishing building a simple connected graph')
        label = graph.sort_edge()  # sort edge in ascending order by the distance
        if not label:
            print('2nd sorting edge did not work, please update code')
            return True
        print('finishing building a simple connected graph')
        Nt = 0
        for i in range(graph.Nnode):
            graph.cid[i] = i  # Reset the cid
        for i in range(graph.Ne):
            v1 = graph.edge[i].id1
            v2 = graph.edge[i].id2
            graph.edge[i].mst_label = False
            graph.edge[i].local_label = False
            graph.edge[i].eid = i

            if graph.cid[v1] == graph.cid[v2]:
                continue
            tree[Nt] = graph.edge[i]
            tmp_cid = graph.cid[v2]
            graph.edge[i].mst_label = True#Used in MST
            Nt += 1
            if MaxCid < tmp_cid:
                MaxCid = tmp_cid
            for j in range(Nt):
                if graph.cid[tree[j].id1] == tmp_cid:
                    graph.cid[tree[j].id1] = graph.cid[v1]
                if graph.cid[tree[j].id2] == tmp_cid:
                    graph.cid[tree[j].id2] = graph.cid[v1]
        graph.Nt = Nt
        graph.tree = tree
        print('cleaning tree finished')
        print('after cleaning Nt=%d, Ne=%d' % (graph.Nt, graph.Ne))
        return graph.cid
        # Clean isolated dens data from mrc


    def build_local_mst(self,graph,point, d2cut,save_path):
        # local MST
        print('local MST building')
        edge_path = os.path.join(save_path, 'edge.txt')
        if os.path.exists(edge_path):
            tmp_edge_data=np.loadtxt(edge_path)
            if len(tmp_edge_data) == graph.Ne:
                for ii in range(graph.Ne):
                    graph.edge[ii].local_label = tmp_edge_data[ii, 6]
            else:
                graph.calcu_local_label(point, d2cut)
        else:
            graph.calcu_local_label(point, d2cut)
        # End local MST
    def setup_tree(self,graph):
        self.Nnode = graph.Nnode
        self.Ne = graph.Nt
        self.Etotal = graph.Ne
        print('after finishing building graph, we get Etotal=%d, Ne=%d' % (self.Etotal, self.Ne))
        if self.Inittree(True):
            return True
        for i in range(graph.Ne):
            if graph.edge[i].mst_label == False:
                continue
            id = graph.edge[i].id1
            eid = graph.edge[i].eid
            self.len += graph.edge[i].d
            self.St = id
            self.ActE[eid] = 1  # set edge's active label
            self.Lpath += 1
        print('MST len:%f' % self.len)
        return False





