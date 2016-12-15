from remove_sm import *
from read_from_file import *
import numpy as npd
import sys
"""
csi_matrix.py

Read Csi matrix and apply some corrections

Naoufal Mahfoudi (c) 2016 mohamed-naoufal.mahfoudi@inria.fr

"""
triangle = [1, 3, 6]
broken_perm = 0


def csi_matrix(csi):
    filename = csi
    c = read_from_file(filename)
    csio = npd.zeros((30,3,2), dtype=complex)
    de_csi = npd.transpose(c[0].csi)
    [m, n, s] = de_csi.shape
    num = npd.size(c)
    op = np.zeros((m, n, s, num), dtype=np.complex_)
    for i in range(0, npd.size(c) - 1):
        perm = c[i].perm
        perm= np.array(perm)
        Nrx = c[i].Nrx
        if Nrx == 1:
            continue
        if np.sum(perm) != triangle[Nrx-1]:
            broken_perm = 1
            print('Invalid CSI File %s with Nrx=%d ',filename,Nrx)
        else:
            c[i].csi=np.array(c[i].csi)
            csio[:, :, :]=(c[i].csi)[:,perm-1,:]
            ret = remove_sm(npd.transpose(csio), c[i].rate)
            op[:, :, :, i] = ret
    return op
