from csi_matrix import *

def phase_correction(csi_corr, csi_target):
    csi_cr = csi_matrix(csi_corr)
    csi_tr = csi_matrix(csi_target)

    rl_dod1 = np.angle(csi_cr[0, 0, 0, 0] * np.conjugate(csi_cr[1, 0, 0, 0]))
    rl_dod2 = np.angle(csi_cr[0, 1, 0, 0] * np.conjugate(csi_cr[1, 1, 0, 0]))
    rl_dod3 = np.angle(csi_cr[0, 2, 0, 0] * np.conjugate(csi_cr[1, 2, 0, 0]))

    rl_doa2 = np.angle(csi_cr[0, 0, 0, 0] * np.conjugate(csi_cr[0, 1, 0, 0]))
    rl_doa3 = np.angle(csi_cr[0, 0, 0, 0] * np.conjugate(csi_cr[0, 2, 0, 0]))

    csi_tro=np.zeros((2, 3,1,np.size(csi_tr, 3)), dtype=complex)
    for k in range(0, np.size(csi_tr, 3) - 1):

        csi_tro[0, 0, 0, k] = csi_tr[0, 0, 0, k]
        csi_tro[0, 1, 0, k] = csi_tr[0, 1, 0, k]
        csi_tro[0, 2, 0, k] = csi_tr[0, 2, 0, k]

        csi_tro[1, 0, 0, k] = csi_tr[1, 0, 0, k] * np.exp(1j * rl_dod1)
        csi_tro[1, 1, 0, k] = csi_tr[1, 1, 0, k] * np.exp(1j * rl_dod2)
        csi_tro[1, 2, 0, k] = csi_tr[1, 2, 0, k] * np.exp(1j * rl_dod3)

        csi_tro[1, 1, 0, k] = csi_tro[1, 1, 0, k] * np.exp(1j * rl_doa2)
        csi_tro[1, 2, 0, k] = csi_tro[1, 2, 0, k] * np.exp(1j * rl_doa3)

        csi_tro[0, 1, 0, k] = csi_tro[0, 1, 0, k] * np.exp(1j * rl_doa2)
        csi_tro[0, 2, 0, k] = csi_tro[0, 2, 0, k] * np.exp(1j * rl_doa3)

    return csi_tro
