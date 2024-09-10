import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

qm = cf.QuantumMechanics(3,xlim=[-15,15],ylim=[-15,15],zlim=[-15,15],
                            xRes=50,yRes=50,zRes=50)
# qm = cf.QuantumMechanics(3,xlim=[-3,3],ylim=[-3,3],zlim=[-3,3],
#                     xRes=50,yRes=50,zRes=50)
qm.psi = qm.calc_hydrogen_state(3,1,1)
qm.psi_f = sp.fft.fftn(qm.psi)
# print(np.real(psi))
# qm.plot_field_in_plane(np.real(psi))
# qm.plot_complex_field_in_plane(psi, normal_vector=[1,0,0])

for n in range(100):
    qm.evolve_schrodinger(10)
    qm.plot_complex_field(qm.psi,phase_blob_threshold=0.3)
    cf.tool_save_plot(n)
cf.tool_make_animation_gif(n)
# qm.plot_complex_field(psi,phase_blob_threshold=0.3)
# plt.show()