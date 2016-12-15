from csi_matrix import *
from numpy import linalg as ln
from numpy import *
from detect_peaks import *
from phase_correction import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
from scipy import ndimage
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
import numpy as np
from matplotlib import cm

def music(csi_corr, csi_target, Ntx, Nrx, d_tx, d_rx,t):

    # csi_corr='csi.dat'
    In = 0
    s = phase_correction(csi_corr, csi_target)

    # sq=np.asarray(squeeze(s[:, :, 0, t:t+2]))

    s_lin =(s[:,:,0,t:t+2].reshape(6,2,order='F'))

    R_hat = np.cov(s_lin)

    D, Q = ln.eig(R_hat)
    '''Sort the eigenvalues in D'''
    Do=np.abs(D)
    D = np.sort(Do)[::-1]
    I = np.argsort(Do)[::-1]
    Q = Q[:, I]

    ''' Compute the Number of signal that are significative'''

    T = np.cumsum(np.real(D))
    for i in range(1, 1, np.size(T)):
        if T(i) >= 0.99 * T(np.size(T)):
            In = i
            break
            #   In = 1
    ''' Get the signal eigenvectors'''
    In=0
    Qs = Q[:,:In]

    ''' Get the noise eigenvectors'''
    Qn = Q[:,In+1:]

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
            s=dot(K.T,Qn)
            music_spectrum[j, k] = 1 / dot(abs(s),abs(s).T)

    ''' compute the mesh and plot the surf of the pseudospectrum '''


    fig = plt.figure()
    ax = fig.gca(projection='3d')
    x = angles2
    y = angles1
    X, Y = np.meshgrid(x, y)
    Z = np.abs(np.squeeze(music_spectrum))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_ylabel('AoA')
    ax.set_xlabel('AoD')
    ax.set_xlim3d(-90, 90)
    ax.set_ylim3d(-90, 90)
    ax.plot_surface(X,Y,Z,rstride=2, cstride=2,cmap=cm.jet,alpha=0.7,linewidth=0.25)


    ''' detect the peaks corresponding to DoD and DoA '''
    detect = detect_peaks(Z)
    index_max = np.column_stack(np.where(detect))
    x_ind = index_max[:, 0]
    y_ind = index_max[:, 1]

    tab = (np.transpose(np.array((Z[x_ind, y_ind], x[x_ind], y[y_ind])))).tolist()
    tab.sort(key=lambda e: e[0], reverse=True)
    myarray = np.asarray(tab[0])
    angles = myarray[1:]
    plt.show()
    return angles