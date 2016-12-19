from csi_matrix import *


def phase_correction(csi_corr, csi_target):
    """phase_correction.py

    Simple Phase correction for one sub-carrier.(Used for the sake of this example.)

    Naoufal Mahfoudi (c) 2016 mohamed-naoufal.mahfoudi@inria.fr

    """

    csi_cr = csi_matrix(csi_corr)
    csi_tr = csi_matrix(csi_target)

    rl_dod = np.zeros(3)
    rl_doa = np.zeros(3)

    for i in range(0, 3):
        rl_dod[i] = np.angle(csi_cr[0, i, 0, 0] * np.conjugate(csi_cr[1, i, 0, 0]))
    for i in range(1, 3):
        rl_doa[i] = np.angle(csi_cr[0, 0, 0, 0] * np.conjugate(csi_cr[0, i, 0, 0]))
    csi_tro = np.zeros((2, 3, 1, np.size(csi_tr, 3)), dtype=complex)
    for k in range(0, np.size(csi_tr, 3) - 1):

        csi_tro[0, :, 0, k] = csi_tr[0, :, 0, k]

        for i in range(0, 3):
            csi_tro[1, i, 0, k] = csi_tr[1, i, 0, k] * np.exp(1j * rl_dod[i])

        for i in range(1, 3):
            csi_tro[:, i, 0, k] = csi_tro[:, i, 0, k] * np.exp(1j * rl_doa[i])
    return csi_tro
