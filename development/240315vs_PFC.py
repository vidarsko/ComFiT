import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import matplotlib.pyplot as plt
import numpy as np

pfc = cf.PhaseFieldCrystal2DSquare(10,10)
# pfc.conf_PFC_from_amplitudes()

pfc.q = 0.95*pfc.q
pfc.conf_PFC_from_amplitudes()
pfc.psi = pfc.psi+ 0.00*np.random.rand(pfc.psi.shape[0],pfc.psi.shape[1])
pfc.psi_f = np.fft.fft2(pfc.psi)

for n in range(50):
    pfc.plot_field(pfc.psi)
    plt.pause(0.01)
    pfc.evolve_PFC(100)
    