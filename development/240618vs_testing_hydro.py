import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import plotly.graph_objs as go
from IPython.display import display, clear_output
import time
import copy


nx = 20
pfc = cf.PhaseFieldCrystal2DTriangular(nx, round(nx/np.sqrt(3)))
eta = pfc.calc_amplitudes_with_dislocation_dipole()
pfc.conf_PFC_from_amplitudes(eta)
pfcH = copy.deepcopy(pfc)



fig = plt.figure(figsize=(15,5))
axs = fig.subplots(2,3)

time_steps = 100

time = []
distance = []
distanceH = []

for n in range(20):

    pfc.evolve_PFC(time_steps)
    time.append(n)
    pfc.plot_field(pfc.psi, ax=axs[0,0], colorbar=False)
    Dnodes = pfc.calc_dislocation_nodes()
    distance.append(np.linalg.norm(np.array(Dnodes[0]['position'])-np.array(Dnodes[1]['position'])))

    pfcH.evolve_PFC_hydrodynamic(time_steps)
    pfcH.plot_field(pfcH.psi[0], ax=axs[1,0], colorbar=False)
    pfcH.plot_vector_field(pfcH.psi[1:2], ax=axs[1,1])
    pfcH.plot_field(np.sqrt(pfcH.psi[1]**2 + pfcH.psi[1]**2), ax=axs[1,2], colorbar=True if n==0 else False)
    DnodesH = pfcH.calc_dislocation_nodes()
    distanceH.append(np.linalg.norm(np.array(DnodesH[0]['position'])-np.array(DnodesH[1]['position'])))

    axs[0,1].plot(time, distance, 'r')
    axs[0,1].plot(time, distanceH, 'b')
    axs[0,1].set_title('Dislocation distance')
    axs[0,1].set_xlabel('Time')
    axs[0,1].set_ylabel('Distance')


    plt.pause(0.01)
    axs[0,0].cla()
    axs[0,1].cla()
    axs[1,0].cla()
    axs[1,1].cla()
    axs[1,2].cla()
