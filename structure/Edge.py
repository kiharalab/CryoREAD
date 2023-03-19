class Edge(object):
    def __init__(self):
        self.d=0
        self.dens=0#Density
        self.id1=0#Node id
        self.id2=0
        self.eid=0#edge id
        self.mst_label=False#the edge is in mst or not
        self.local_label=False#Whether it is in local region
        #if there exists a node is in acceptable range to an edge's two nodes, then the edge will be marked as local_label=True
        self.keep_label=False#the 100% will be on the tree will be marked as keep_label=True
    def copy_edge(self,in_edge):
        self.d=in_edge.d
        self.dens=in_edge.dens
        self.id1=in_edge.id1
        self.id2=in_edge.id2
        self.eid=in_edge.eid
        self.mst_label=in_edge.mst_label
        self.local_label=in_edge.local_label
        self.keep_label=in_edge.keep_label

