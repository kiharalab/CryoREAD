

import numpy as np

def calculate_distance(x,y):
    distance=0
    for k in range(3):
        distance+=(x[k]-y[k])**2
    distance = np.sqrt(distance)
    return distance

def calculate_cosine_value(edge1,edge2,edge3):
    nominator = edge1**2+edge2**2-edge3**2
    denominator = 2*edge1*edge2
    return nominator/denominator
