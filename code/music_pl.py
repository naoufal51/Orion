"""
music_pl.py
Joint AoD AoA estimation using 2D music Algorithm
Naoufal Mahfoudi (c) 2016 mohamed-naoufal.mahfoudi@inria.fr
"""
from numpy import linalg as ln
from numpy import *
from detect_peaks import *
from phase_correction import *
import numpy as np


def music_pl(s, Ntx, Nrx, d_tx, d_rx, t):

    """Joint AoD AoA estimation using 2D music Algorithm
        using raw.dat file as Input

      Args:
        s: Phase corrected csi data
        Ntx: Number of transmitting antennas.
        Nrx: Number of receiving antennas.
        d_tx: distance between the Tx antennas.
        d_rx: distance between the Rx antennas.
        t: The packet index number

      Returns:
        An array 'angles' with the estimated DoA and DoD in this order.

      """
    s_lin = (s[:, :, 0, t:t + 2].reshape(6, 2, order='F'))

    '''Compute the covariance matrix and the eigendecompositon'''
    R_hat = np.cov(s_lin)
    D, Q = ln.eig(R_hat)

    '''Sort the eigenvalues in D'''
    Do = np.abs(D)
    D = np.sort(Do)[::-1]
    I = np.argsort(Do)[::-1]
    Q = Q[:, I]

    ''' Compute the Number of signal that are significative'''
    T = np.cumsum(np.real(D))
    for i in range(1, 1, np.size(T)):
        if T(i) >= 0.99 * T(np.size(T)):
            In = i
            break

    ''' Get the signal eigenvectors'''
    In = 0
    Qs = Q[:, :In]

    ''' Get the noise eigenvectors'''
    Qn = Q[:, In + 1:]

    ''' Angles at which MUSIC Pseudospectrum  will be computed '''
    angles1 = np.arange(-90, 90, 1)
    angles2 = np.arange(-90, 90, 1)

    '''Compute steering vectors corresponding values in angles'''
    a1 = np.exp(-1.j * 2 * np.pi * d_rx * np.tensordot(arange(Nrx), sin(angles1 * np.pi / 180), 0))
    a2 = np.exp(-1.j * 2 * np.pi * d_tx * np.tensordot(arange(Ntx), sin(angles1 * np.pi / 180), 0))

    '''Compute MUSIC "spectrum" '''
    music_spectrum = np.zeros((np.size(angles1), np.size(angles2)), dtype=complex)
    for k in range(1, np.size(angles2)):
        for j in range(1, np.size(angles1)):
            K = np.kron(a1[:, j], a2[:, k])
            s = dot(K.T, Qn)
            music_spectrum[j, k] = 1 / dot(abs(s), abs(s).T)

    ''' compute the mesh and plot the surf of the pseudospectrum '''
    x = angles2
    y = angles1
    X, Y = np.meshgrid(x, y)
    Z = np.abs(np.squeeze(music_spectrum))

    ''' detect the peaks corresponding to DoD and DoA '''
    detect = detect_peaks(Z)
    index_max = np.column_stack(np.where(detect))
    x_ind = index_max[:, 0]
    y_ind = index_max[:, 1]
    tab = (np.transpose(np.array((Z[x_ind, y_ind], x[x_ind], y[y_ind])))).tolist()
    tab.sort(key=lambda e: e[0], reverse=True)
    myarray = np.asarray(tab[0])
    angles = myarray[1:]

    return angles
