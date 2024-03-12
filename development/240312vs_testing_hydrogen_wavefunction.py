import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt

qm = cf.QuantumMechanics(3,xlim=[-10,10],ylim=[-10,10],zlim=[-10,10],
                            xRes=150,yRes=150,zRes=150)
psi = qm.calc_hydrogen_state(2,1,1)

print(psi)
qm.plot_field_in_plane(np.real(psi))
plt.show()