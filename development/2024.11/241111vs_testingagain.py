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


nx=10
B_x=0.98
psi0=0
t=1/2/B_x
v=1/3/B_x
Delta_B = -0.14
time_limit = 100000
time_step = 1000

pfc0 = cf.PhaseFieldCrystal2DTriangular(nx,round(nx/np.sqrt(3)), t=t, r=Delta_B/B_x, v=v, psi0=psi0)
pfc0.conf_PFC_from_amplitudes()  

res = 101
strains = np.linspace(-0.2,0.2,21)
time_to_nucleation = np.nan*np.zeros_like(strains)


strain_counter = 0
finished = False
from_low_strain=True
for m in range(len(strains)):
    if not finished:
        strain = strains[strain_counter] if from_low_strain else strains[-strain_counter-1]
        j = strain_counter if from_low_strain else -strain_counter-1
        strain_counter += 1
        print("Considering strain = ", strains[j])
        
        pfc = copy.deepcopy(pfc0)
        distortion = np.eye(2)*strain
        pfc.conf_apply_distortion(distortion, update_q_and_a_vectors=False)
        pfc.psi = pfc.psi + 1e-2*(np.random.rand(pfc.psi.shape[0], pfc.psi.shape[1])-0.5)
        pfc.psi_f = sp.fft.fftn(pfc.psi)

        free_energy  = pfc.calc_free_energy()

        free_energy_sharp_drop=False
        free_energy_stable=False

        time = 0
        free_energy = pfc.calc_free_energy()
        initial_free_energy = free_energy
        while not free_energy_stable and not free_energy_sharp_drop:
            pfc.evolve_PFC(100,suppress_output=True)
            time += 10
            free_energy_new = pfc.calc_free_energy()

            if (free_energy - free_energy_new) < 1e-14:
                free_energy_stable=True
                print("Free energy stable")
                if from_low_strain:
                    from_low_strain=False
                    strain_counter=0
                else:
                    finished = True
                    
            
            if free_energy_new/initial_free_energy < 0.9:
                free_energy_sharp_drop=True
                print("Free energy sharp drop")
            
            free_energy = free_energy_new

        time_to_nucleation[j] = time
    
    

        
plt.plot(strains, time_to_nucleation)
plt.show()
    