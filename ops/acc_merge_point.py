
from numba import jit
import numpy as np
import math
@jit(nogil=True,nopython=True)
def acc_merge_point(Ncd,dens,dmin,rv_range,rdcut,stock,cd,d2cut,member):
    if True:
        for i in range(Ncd-1):
            if i%10000==0:
                print(i)
            tmp=np.zeros(3)
            if (dens[i]-dmin)*rv_range < rdcut:
                stock[i]=0#Label the small density parts as unused parts
            if stock[i]==0:
                continue
            for j in range(i+1,Ncd):
                if stock[j]==0:
                    continue
                d2=0
                for k in range(3):
                    tmp[k]=cd[i][k]-cd[j][k]
                    d2+=tmp[k]**2
                if d2<d2cut:
                    #Mark the merged points to where it goes
                    if dens[i]>dens[j]:
                        stock[j]=0
                        member[j]=i
                    else:
                        stock[i]=0
                        member[i]=j
                        break#jump out of the second rotation, since i has been merged
        #Update member data, to updata some son/grandson points to original father point
        for i in range(Ncd):
            now=int(member[i])
            while now!=member[now]:#If it's not merged points, it will totates to find the father point(merged point)
                now=int(member[now])
            member[i]=now
        return stock,member
