import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp


pfc = cf.PhaseFieldCrystal2DTriangular(30, 18)

# This creates a standard orientation of the crystal
pfc.conf_PFC_from_amplitudes()

# Create the rotated field
psi_rotated = pfc.calc_PFC_from_amplitudes(rotation=np.pi/6)

# Specify the region centered at the mid position with radius 5 a0.
inclusion_region = pfc.calc_region_disk(pfc.rmid, 10*pfc.a0)
pfc.evolve_PFC(100) #Smooth the interface

# Set the rotated field in the inclusion region
pfc.psi[inclusion_region] = psi_rotated[inclusion_region]
pfc.psi_f = sp.fft.fftn(pfc.psi)

for n in range(100):
    pfc.evolve_PFC_hydrodynamic(100)
    pfc.plot_field(pfc.psi[0])
    cf.tool_save_plot(n)
cf.tool_make_animation_gif(n,'Inclusion_PFC_hydrodynamic')