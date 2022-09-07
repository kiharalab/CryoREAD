import numpy as np
class Node(object):
    def __init__(self):
        self.N=0#number of the connected edges of the node
        self.e=[]#edge list
        self.cid=0#chain id of the node
        self.node_id=0#
        #["O", "C", "N", "P", "others"]
        self.atom_prob = np.zeros(5)
        self.log_atom = np.zeros(5)
        self.count_atom = 0
        # ["sugar", "phosphorus", "A", "UT", "C", "G", "protein"]
        self.nuc_prob = np.zeros(7)
        self.count_nuc = 0
        self.log_nuc = np.zeros(7)
        # main suggest "sugar", "phosphorus"
        # side suggest  "A", "UT", "C", "G",
        #chain_label_list = ['DRNA_main',"DRNA_side","protein"]
        self.chain_prob = np.zeros(3)
        self.count_chain = 0
        self.log_chain = np.zeros(3)
        #base predictions
        self.base_prob = np.zeros(5)
        self.count_base = 0
        self.log_base = np.zeros(5)
