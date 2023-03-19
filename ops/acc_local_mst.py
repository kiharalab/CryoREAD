from numba import jit
import numpy as np
import math

@jit(nogil=True,nopython=True)
def acc_local_mst(cid,Nnode,Ne,id1,id2,merged_cd_dens,d2cut,local_label):
    """
    :param cid:
    :param Nnode:
    :param Ne:
    :param id1:
    :param id2:
    :param merged_cd_dens:
    :param d2cut:
    :param local_label:
    :return:
    purpose: if there exists a node is in acceptable range to an edge's two nodes, then the edge will be marked as local_label=True
    """
    if True:

        count=0
        for i in range(Nnode):
            vec=np.zeros(3)
            #if i%100==0:
            #    print(i)
            for j in range(Nnode):
                cid[j]=j
            for k in range(Ne):
                v1=int(id1[k])
                v2=int(id2[k])
                if cid[v1]==cid[v2]:
                    continue
                #dist
                d2=0
                for l in range(3):
                    vec[l]=merged_cd_dens[i,l]-merged_cd_dens[v1,l]
                    d2+=vec[l]**2
                if d2>d2cut:
                    continue
                d2=0
                for l in range(3):
                    vec[l]=merged_cd_dens[i,l]-merged_cd_dens[v2,l]
                    d2+=vec[l]**2
                if d2>d2cut:
                    continue
                local_label[k]=1
                count+=1
                tmpid=cid[v2]
                for l in range(Nnode):
                    if cid[l]==tmpid:#Update cid
                        cid[l]=cid[v1]
    #print(count)
    return local_label,count
