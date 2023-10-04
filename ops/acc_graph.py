from numba import jit
import numpy as np
import math

@jit(nogil=True,nopython=True)
def acc_new_adj_matrix(Nori,merged_data,adj,origrid,Nnode):
    Ne=0
    abanden=0
    check_set = set()
    for ii in range(Nori):
        m1 = int(merged_data[ii, 0])
        merged_id1 = int(merged_data[m1, 5])
        if ii % 10000 == 0:
            print(ii)
        if merged_id1 == -1:
            continue
        for jj in range(ii + 1, Nori):
            m2 = int(merged_data[jj, 0])
            if m1 == m2:
                continue

            merged_id2 = int(merged_data[m2, 5])
            if merged_id2 == -1 or merged_id1 == merged_id2:  # To avoid that point with stock=0, very lower density, however, never be combined to others
                continue
            if merged_id1*Nnode+merged_id2 in check_set:
                continue
            adjunct_label = False
            for kk in range(3):
                if (origrid[ii][kk] - origrid[jj][kk]) ** 2 > 1:
                    adjunct_label = True  # It will not connected anymore
                    break
            if adjunct_label:
                abanden += 1
                continue
            check_set.add(merged_id1 * Nnode + merged_id2)
            adj[Ne] = merged_id1*Nnode+merged_id2

            Ne += 1  # count for effective edges
    return adj,Ne
@jit(nogil=True,nopython=True)
def acc_adj_matrix(Nori,merged_data,adj,origrid):
    if True:
        Ne=0
        abanden=0
        #count_useful=0
        for ii in range(Nori):#original point before merging
            m1=int(merged_data[ii,0])
            merged_id1=int(merged_data[m1,5])
            if ii%10000==0:
                print(ii)
            if merged_id1==-1:
                continue
            for jj in range(ii+1,Nori):
                m2=int(merged_data[jj,0])
                if m1==m2:
                    continue

                merged_id2=int(merged_data[m2,5])
                if adj[merged_id1,merged_id2]:
                    continue
                if merged_id2==-1 or merged_id1==merged_id2:#To avoid that point with stock=0, very lower density, however, never be combined to others
                    continue

                adjunct_label=False
                for kk in range(3):
                    if (origrid[ii][kk]-origrid[jj][kk])**2>1:
                        adjunct_label=True#It will not connected anymore
                        break
                if adjunct_label:
                    abanden+=1
                    continue
                adj[merged_id1][merged_id2]=1
                adj[merged_id2][merged_id1]=1
                Ne+=1#count for effective edges

    return adj,Ne

@jit(nogil=True,nopython=True)
def acc_edge_prob(Ne,id1,id2,merged_cd_dens,dens,edge_d,fsiv,fmaxd,xdim,ydim,zdim,mrc_dense):
    if True:
        for ii in range(Ne):
            if ii%1000==0:
                print(ii)
            v1=int(id1[ii])
            v2=int(id2[ii])
            vec=np.zeros(3)
            for kk in range(3):
                vec[kk]=merged_cd_dens[v2,kk]-merged_cd_dens[v1,kk]
            MinDens=merged_cd_dens[v1,3]
            cd1=np.zeros(3)
            for jj in range(11):
                for kk in range(3):
                    cd1[kk]=merged_cd_dens[v1,kk]+vec[kk]*0.1*jj
                tmp_dens=mean_shift_pos(cd1,fsiv,fmaxd,xdim,ydim,zdim,mrc_dense)
                if tmp_dens<MinDens:
                    MinDens=tmp_dens
            if MinDens<=1e-4:#very rare case when it come to border of the box
                MinDens=1e-4 #to avoid very big values #ignore those boundary cases#(merged_cd_dens[v1,3]+merged_cd_dens[v2,3])/2
            dens[ii]=edge_d[ii]/MinDens#prob bigger, distance smaller is our favorite choice
    return dens

@jit(nogil=True,nopython=True)
def acc_edge_density(Ne,id1,id2,merged_cd_dens,dens,edge_d,fsiv,fmaxd,xdim,ydim,zdim,mrc_dense):
    if True:
        for ii in range(Ne):
            if ii%1000==0:
                print(ii)
            v1=int(id1[ii])
            v2=int(id2[ii])
            vec=np.zeros(3)
            for kk in range(3):
                vec[kk]=merged_cd_dens[v2,kk]-merged_cd_dens[v1,kk]
            MinDens=merged_cd_dens[v1,3]#Notice 09/25/2018, important update, we need to get the unnormalized density here, actually we did
            cd1=np.zeros(3)
            for jj in range(11):
                for kk in range(3):
                    cd1[kk]=merged_cd_dens[v1,kk]+vec[kk]*0.1*jj
                tmp_dens=mean_shift_pos(cd1,fsiv,fmaxd,xdim,ydim,zdim,mrc_dense)
                if tmp_dens<MinDens:
                    MinDens=tmp_dens
            dens[ii]=MinDens*edge_d[ii]
    return dens
@jit(nogil=True,nopython=True)
def mean_shift_pos(pos,fsiv,fmaxd,xdim,ydim,zdim,dens):
        #set up filter

        #mean shifting
        stp=np.zeros(3)
        endp=np.zeros(3)
        #start point
        for k in range(3):
            stp[k]=(int)(pos[k]-fmaxd)
            if stp[k]<0:
                stp[k]=0
            endp[k]=(int)(pos[k]+fmaxd+1)
        #print(endp[0])
        #print(xdim)
        if endp[0]>=xdim:
            endp[0]=xdim
        if endp[1]>=ydim:
            endp[1]=ydim
        if endp[2]>=zdim:
            endp[2]=zdim
        dtotal=0
        for xp in range(int(stp[0]),int(endp[0])):
            rx=(float)(xp-pos[0])**2
            for yp in range(int(stp[1]),int(endp[1])):
                ry=(float)(yp-pos[1])**2
                for zp in range(int(stp[2]),int(endp[2])):
                    rz=(float)(zp-pos[2])**2
                    d2=rx+ry+rz
                    v=math.exp(-1.5*d2*fsiv)*dens[xp,yp,zp]
                    dtotal+=v
        return dtotal
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

@jit(nogil=True,nopython=True)
def  acc_get_density_cut_twographextra(check_list,merged_data,adj,density_record,origrid,Nnode,
                                 cutoff_index,origrid_cutoff_index):
    """
    :param Nori: number of grid points that can be used
    :param merged_data: an array with format # init_id, x,y,z,density, merged_to_id
    :param adj: used to save the adjacency info
    :param density_record: used to save the edge density information
    :param origrid: origin grid array
    :param Nnode: number of merged node reamined
    :param cutoff_index: boundary of two nodes
    :return:
    #only focused pho or sugar internal connections
    """
    Ne = 0
    abanden = 0
    for ii in check_list:
        ii = int(ii)
        m1 = int(merged_data[ii, 0])#original id
        if ii<origrid_cutoff_index:
            merged_id1 = int(merged_data[m1, 5])#link to  merged point id
        else:
            merged_id1 = int(merged_data[m1+origrid_cutoff_index, 5])
        if ii>=origrid_cutoff_index:
            merged_id1 +=cutoff_index
        if ii % 10000 == 0:
            print(ii)
        if merged_id1 == -1:#Merged node, ignore
            continue
        for jj in check_list:
            jj = int(jj)
            if ii==jj:
                continue
            m2 = int(merged_data[jj, 0])
            if ii<origrid_cutoff_index and jj>=origrid_cutoff_index:
                continue
            if ii>=origrid_cutoff_index and jj<origrid_cutoff_index:
                continue
            if jj>=origrid_cutoff_index:
                merged_id2 = int(merged_data[m2+origrid_cutoff_index, 5])
            else:
                merged_id2 = int(merged_data[m2, 5])
            if merged_id2 == -1:
                continue

            if jj>=origrid_cutoff_index:
                merged_id2+=cutoff_index#map to real merge id
            if merged_id1==merged_id2:
                continue
            adjunct_label = False
            for kk in range(3):
                if (origrid[ii][kk] - origrid[jj][kk]) ** 2 > 1:
                    adjunct_label = True  # It will not connected anymore
                    break
            if adjunct_label:
                abanden += 1
                continue
            adj[Ne] = merged_id1 * Nnode + merged_id2#Redundant, every possible connect save twice
            density_record[Ne]=min(merged_data[ii, 4],merged_data[jj, 4])#Once one don't exist, then the connection will be not exist by this specific ori points
            Ne += 1  # count for effective edges
    return adj, Ne,density_record




@jit(nogil=True,nopython=True)
def acc_get_density_cut_twograph(Nori,merged_data,adj,density_record,origrid,Nnode,
                                 cutoff_index,origrid_cutoff_index):
    """
    :param Nori: number of grid points that can be used
    :param merged_data: an array with format # init_id, x,y,z,density, merged_to_id
    :param adj: used to save the adjacency info
    :param density_record: used to save the edge density information
    :param origrid: origin grid array
    :param Nnode: number of merged node reamined
    :param cutoff_index: boundary of two nodes
    :return:
    """
    Ne = 0
    abanden = 0
    for ii in range(origrid_cutoff_index):
        m1 = int(merged_data[ii, 0])#original id
        merged_id1 = int(merged_data[m1, 5])#link to  merged point id

        if ii % 10000 == 0:
            print(ii)
        if merged_id1 == -1:#Merged node, ignore
            continue
        for jj in range(origrid_cutoff_index, Nori):
            m2 = int(merged_data[jj, 0])
            merged_id2 = int(merged_data[m2+origrid_cutoff_index, 5])
            if merged_id2 == -1:
                continue
            merged_id2+=cutoff_index#map to real merge id
            assert merged_id1 != merged_id2#it should be impossible to be equal
            adjunct_label = False
            for kk in range(3):
                if (origrid[ii][kk] - origrid[jj][kk]) ** 2 > 1:
                    adjunct_label = True  # It will not connected anymore
                    break
            if adjunct_label:
                abanden += 1
                continue
            adj[Ne] = merged_id1 * Nnode + merged_id2#Redundant, every possible connect save twice
            density_record[Ne]=min(merged_data[ii, 4],merged_data[jj, 4])#Once one don't exist, then the connection will be not exist by this specific ori points
            Ne += 1  # count for effective edges
    return adj, Ne,density_record

@jit(nogil=True,nopython=True)
def acc_get_density_cut(Nori,merged_data,adj,density_record,origrid,Nnode):
    """
    :param Nori: number of grid points that can be used
    :param merged_data: an array with format # init_id, x,y,z,density, merged_to_id
    :param adj: used to save the adjacency info
    :param density_record: used to save the edge density information
    :param origrid: origin grid array
    :param Nnode: number of merged node reamined
    :return:
    """
    Ne = 0
    abanden = 0
    for ii in range(Nori):
        m1 = int(merged_data[ii, 0])#original id
        merged_id1 = int(merged_data[m1, 5])#link to  merged point id
        if ii % 10000 == 0:
            print(ii)
        if merged_id1 == -1:#Merged node, ignore
            continue
        for jj in range(ii + 1, Nori):
            m2 = int(merged_data[jj, 0])
            if m1 == m2:
                continue

            merged_id2 = int(merged_data[m2, 5])
            if merged_id2 == -1 or merged_id1 == merged_id2:  # To avoid that point with stock=0, very lower density, however, never be combined to others
                continue

            adjunct_label = False
            for kk in range(3):
                if (origrid[ii][kk] - origrid[jj][kk]) ** 2 > 1:
                    adjunct_label = True  # It will not connected anymore
                    break
            if adjunct_label:
                abanden += 1
                continue
            #make sure init points has connections
            adj[Ne] = merged_id1 * Nnode + merged_id2#Redundant, every possible connect save twice
            density_record[Ne]=min(merged_data[ii, 4],merged_data[jj, 4])#Once one don't exist, then the connection will be not exist by this specific ori points
            Ne += 1  # count for effective edges
    return adj, Ne,density_record
