from sm_matrices import *

def remove_sm(csi, rate):
    m = np.size(csi, 0)
    n = np.size(csi, 1)
    s = np.size(csi, 2)
    ret = np.zeros((m, n, s), dtype=np.complex_)

    if m == 1:
        ret = csi
        return ret
    sm = []
    cond=(np.bitwise_and(rate, 2048) == 2048)

    if cond:
        if m == 3:
            sm = sm_3_40
        elif m == 2:
            sm = sm_2_40
    else:
        if m == 3:
            sm = sm_3_20
        elif m == 2:
            sm = sm_2_20
    for i in range(0, s):
        t = np.array(csi)[:, :, i]
        ret[:, :, i] = np.transpose(np.transpose(t) * np.transpose(sm))
    return ret
#         ret[:, :, i] = np.matrix(np.squeeze(t)).T * np.matrix(sm_2_20).T
