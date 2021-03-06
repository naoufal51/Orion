"""kalman_fil.py

Classic Kalman filter for statistical noise mitigation.

Naoufal Mahfoudi (c) 2016 mohamed-naoufal.mahfoudi@inria.fr

"""
import pandas as pd

from hampel import *


def kalman_fil(csv):
    """Kalman filter

      Args:
        csv: csv file containing the AoA and AoD.

      Returns:
        array of filter DoA estimations.
        array of filter DoD estimations.

      """
    df = pd.read_csv(csv, sep=',', header=None)
    angles = array(df.values)
    angles1 = hampel(angles[0, :], m=2)
    angles2 = hampel(angles[1, :], m=2)

    dt = 0.1

    A = matrix([[1, 0, dt, 0],
                [0, 1, 0, dt],
                [0, 0, 1, 0],
                [0, 0, 0, 1]])
    B = matrix([[(dt ** 2 / 2)], [(dt ** 2 / 2)], [dt], [dt]])
    u = .005
    H = matrix([[1, 0, 0, 0], [0, 1, 0, 0]])
    Tz = matrix([[1, 0],
                 [0, 1]])
    Acc_ns = 0.03
    Tx = matrix([[dt ** 4 / 4, 0, dt ** 3 / 2, 0],
                 [0, dt ** 4 / 4, 0, dt ** 3 / 2],
                 [dt ** 3 / 2, 0, dt ** 2, 0],
                 [0, dt ** 3 / 2, 0, dt ** 2]]) * Acc_ns ** 2
    P = Tx

    pos_estimate = matrix([[angles[0, 0]], [angles[1, 0]], [0], [0]])
    pos_loc_meas_p = np.zeros((2, size(angles, 1)))
    DOD = np.zeros(size(angles, 1))
    DOA = np.zeros(size(angles, 1))

    for t in range(1, size(angles, 1)):
        '''Measurement and parameter estimation'''
        pos_loc_meas_p[0, t] = angles1[t]
        pos_loc_meas_p[1, t] = angles2[t]

        pos_loc_meas = pos_loc_meas_p[:, t]
        '''Predict Stage'''
        pos_estimate = A * pos_estimate + B * u
        P = A * P * A.T + Tx

        '''Update Stage'''
        y = pos_loc_meas - H * pos_estimate
        s = matrix(H * P * H.T + Tz)  # residual covariance
        K = P * H.T * s.I
        pos_estimate = pos_estimate + K * y
        I = matrix(eye(A.shape[0]))
        P = (I - K * H) * P

        '''Storing Stage'''
        DOA[t] = (pos_estimate[0, 0])
        DOD[t] = (pos_estimate[0, 1])
    return DOA, DOD, angles1, angles2
