import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp


pfc = cf.PhaseFieldCrystal3DSimpleCubic(3,3,3)

pfc.conf_PFC_from_amplitudes(rotation=[0,0.0,0.0])

orientation_field = pfc.calc_orientation_field()
# pfc.plot_field(pfc.psi)
# pfc.plot_orientation_field(orientation_field)
orientation_field_norm = np.sqrt(orientation_field[0]**2 + orientation_field[1]**2 + orientation_field[2]**2 + orientation_field[3]**2)
theta = np.arccos(orientation_field[0]/orientation_field_norm)
nnorm = np.sqrt(orientation_field[1]**2 + orientation_field[2]**2 + orientation_field[3]**2)
nx = orientation_field[1]/nnorm
ny = orientation_field[2]/nnorm 
nz = orientation_field[3]/nnorm 
print(theta)
vector_field = theta*np.array([nx,ny,nz])

# print(vector_field)
pfc.plot_vector_field(vector_field, spacing=5)

plt.show()