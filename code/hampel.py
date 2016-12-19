"""
hampel.py

Outlier detection.

Naoufal Mahfoudi (c) 2016 mohamed-naoufal.mahfoudi@inria.fr

"""
from music import *

def hampel(data, m=2.):
    """hampel outlier detector

       Args:
         data: 1D array containing the estimates

       Returns:
         outlier mitigated array.
       """
    res = zeros(size(data))
    for t in range(5, size(data)):
        d = np.abs(data[t] - np.median(data[t - 5:t + 5]),dtype=float)
        mad = np.median(np.abs(data[t - 5:t + 5] - np.median(data[t - 5:t + 5])))
        a = d/mad
        if a >0.1:
            res[t] = np.median(data[t - 5:t + 5])
        else:
            res[t] = data[t]
    return res
