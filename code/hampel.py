"""
hampel.py

Outlier detection.

Naoufal Mahfoudi (c) 2016 mohamed-naoufal.mahfoudi@inria.fr

"""
from music import *

def hampel(data, m=2.):
    res = zeros(size(data))
    for t in range(5, size(data)):
        d = np.abs(data[t] - np.median(data[t - 3:t + 3]))
        mad = np.median(np.abs(data[t - 3:t + 3] - np.median(data[t - 3:t + 3])))
        a = d / mad
        if a > m:
            res[t] = np.median(data[t - 3:t + 3])
        else:
            res[t] = data[t]
    return res