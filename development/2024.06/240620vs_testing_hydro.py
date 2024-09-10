import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import plotly.graph_objs as go
import time
import copy


nx = 30
# pfc = cf.PhaseFieldCrystal2DTriangular(nx, round(nx/np.sqrt(3)))
pfc = cf.PhaseFieldCrystal2DSquare(nx, nx)
eta = pfc.calc_amplitudes_with_dislocation_dipole()
pfc.conf_PFC_from_amplitudes(eta)
pfc.evolve_PFC(100)
pfcH = copy.deepcopy(pfc)




time_steps = 10

time = []
distance = []
distanceH = []

rho0exp = -6
gammaSexp = -6


fig = plt.figure(figsize=(15,5))
axs = fig.subplots(2,3)
for n in range(200):

    pfc.evolve_PFC(time_steps)
    time.append(n)
    pfc.plot_field(pfc.psi, ax=axs[0,0], colorbar=False)
    Dnodes = pfc.calc_dislocation_nodes()
    if len(Dnodes) == 2:
        distance.append(np.linalg.norm(np.array(Dnodes[0]['position'])-np.array(Dnodes[1]['position'])))
    else:
        distance.append(0)

    pfcH.evolve_PFC_hydrodynamic(time_steps, gamma_S = 2**gammaSexp,
                                            rho0 = 2**rho0exp)
    pfcH.plot_field(pfcH.psi[0], ax=axs[1,0], colorbar=False)
    pfcH.plot_vector_field(pfcH.psi[1:2], ax=axs[1,1])
    pfcH.plot_field(np.sqrt(pfcH.psi[1]**2 + pfcH.psi[1]**2), ax=axs[1,2], colorbar=True)
    DnodesH = pfcH.calc_dislocation_nodes()
    if len(DnodesH) == 2:
        distanceH.append(np.linalg.norm(np.array(DnodesH[0]['position'])-np.array(DnodesH[1]['position'])))
    else:
        distanceH.append(0)
    
    axs[0,1].plot(time, distance, 'r')
    axs[0,1].plot(time, distanceH, 'b')
    axs[0,1].set_title('Dislocation distance')
    axs[0,1].set_xlabel('Time')
    axs[0,1].set_ylabel('Distance')



    # plt.pause(0.01)
    cf.tool_save_plot(n,image_size_inches=(15,5))
    
    fig.clear()
    fig = plt.figure(figsize=(15,5))
    axs = fig.subplots(2,3)

cf.tool_make_animation_gif(n, name = f'rho0_2_{rho0exp}_gammaS_2_{gammaSexp}')
