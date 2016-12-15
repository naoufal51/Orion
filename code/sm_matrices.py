"""
sm_matrices.py

sm matrices used by the intel 5300 agn

Naoufal Mahfoudi (c) 2016 mohamed-naoufal.mahfoudi@inria.fr

"""
import numpy as np

sm_1 = 1

sm_2_20 = np.matrix('1 1; 1 -1') / np.sqrt(2)

sm_2_40 = np.matrix('1 1; 1j 1') / np.sqrt(2)

sm_3_20 = np.matrix('-2*np.pi/16 -2*np.pi/(80/33) 2*np.pi/(80/3);'
                    '2*np.pi/(80/23) 2*np.pi/(48/13) 2*np.pi/(240/13);'
                    '-2*np.pi/(80/13) 2*np.pi/(240/37) 2*np.pi/(48/13)')
sm_3_20 = np.exp(1j * sm_3_20) / np.sqrt(3)

sm_3_40 = np.matrix('-2*np.pi/16 -2*np.pi/(80/13) 2*np.pi/(80/23);'
                    '2*np.pi/(80/37) 2*np.pi/(48/11) 2*np.pi/(240/107);'
                    '-2*np.pi/(80/7) 2*np.pi/(240/83) 2*np.pi/(48/11)')
sm_3_20 = np.exp(1j * sm_3_40) / np.sqrt(3)