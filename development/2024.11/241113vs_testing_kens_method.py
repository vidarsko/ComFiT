import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scipy as sp
import copy
import os
from datetime import datetime
import json 
import pickle


B_x=0.98
psi0=0
t=1/2/B_x
v=1/3/B_x
Delta_B = -0.14

res = 101
strains = np.linspace(-0.2,0.2,81)
time_to_nucleation = np.nan*np.zeros_like(strains)

strain_counter = 0
finished = False
from_low_strain=True
for m in range(len(strains)):
    if not finished:
        strain = strains[strain_counter] if from_low_strain else strains[-strain_counter-1]
        j = strain_counter if from_low_strain else -strain_counter-1
        strain_counter += 1
        print(f"Considering strain = {strains[j]:.3f}")

        xmax = 1000*0.63
        ax = 4*np.pi/(np.sqrt(3)*(1-strain))
        ay = 4*np.pi/(1-strain)

        nx = int(xmax/ax + 0.5)
        ny = int(xmax/ay + 0.5)

        dx = nx*ax/1000
        dy = ny*ay/1000

        micro_resolution = ax/dx 

        pfc0 = cf.PhaseFieldCrystal2DTriangular(nx,ny, t=t, r=Delta_B/B_x, v=v, psi0=psi0)
        pfc0.conf_PFC_from_amplitudes()  

        pfc = copy.deepcopy(pfc0)
        distortion = np.eye(2)*strain
        pfc.conf_apply_distortion(distortion, update_q_and_a_vectors=False)
        pfc.psi = pfc.psi + 1e-2*(np.random.rand(pfc.psi.shape[0], pfc.psi.shape[1])-0.5)
        pfc.psi_f = sp.fft.fftn(pfc.psi)
        pfc.evolve_PFC(100,suppress_output=True)

        free_energy  = pfc.calc_free_energy()

        free_energy_sharp_drop=False
        free_energy_stable=False

        time = 0
        free_energy = pfc.calc_free_energy()
        initial_free_energy = free_energy
        while not free_energy_stable and not free_energy_sharp_drop:
            pfc.evolve_PFC(100,suppress_output=True)
            time += 10
            time_to_nucleation[j] = time
            free_energy_new = pfc.calc_free_energy()

            if (free_energy - free_energy_new) < 1e-15:
                free_energy_stable=True
                print("Stable: Free energy stable")
                time_to_nucleation[j] = np.nan
                if from_low_strain:
                    pfc_low_strain_stable = copy.deepcopy(pfc)
                    from_low_strain=False
                    strain_counter=0
                else:
                    pfc_high_strain_stable = copy.deepcopy(pfc)
                    finished = True
            
            if (initial_free_energy-free_energy_new)/abs(initial_free_energy) > 0.1:
                free_energy_sharp_drop=True
                print("Unstable: Free energy sharp drop")
                if from_low_strain:
                    pfc_low_strain_unstable = copy.deepcopy(pfc)
                else:
                    pfc_high_strain_unstable = copy.deepcopy(pfc)
            
            free_energy = free_energy_new

        
    
    

        
plt.plot(strains, time_to_nucleation)
plt.axvline(x=-0.1425, color='r', linestyle='--')
plt.axvline(x=0.1175, color='r', linestyle='--')
plt.show()
    