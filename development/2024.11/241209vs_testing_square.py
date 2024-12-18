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
from tqdm import tqdm

data_file_path = 'C:/Users/vidar/Desktop/ComFiT/development/2024.10/stab7mode_98.dat'
data = np.loadtxt(data_file_path)
# Save column one into strain_ana and column two into deltab_ana
strain_ana = data[:, 0]
deltab_ana = data[:, 1]

nx=20
B_x=0.98
psi0=0
t=1/2/B_x
v=1/3/B_x

time_step=100

Delta_Bs = np.arange(-0.15, 0.07, 0.005)
strains = np.arange(-0.2, 0.205, 0.005)

Delta_Bs_mesh, strains_mesh = np.meshgrid(Delta_Bs,strains)
time_to_nucleation = np.nan*np.zeros_like(Delta_Bs_mesh)

def calc_derivatives(F, time_step):
    n = len(F)
    derivatives = np.zeros(n-1)
    
    for i in range(n-1):
        dFdt = np.gradient(F, time_step)
        derivatives[i] = dFdt[0]

        F = np.zeros(n-i)
        F = dFdt[:-1]

    return derivatives



def update_plot():
    plt.clf()
    plt.pcolormesh(strains_mesh.transpose(), Delta_Bs_mesh.transpose(), np.log10(time_to_nucleation.transpose()), shading='auto')
    plt.colorbar(label='Time to Nucleation (powers of 10)')
    plt.plot(strain_ana, deltab_ana, 'r', label='Analytical')
    plt.xlabel('Strain')
    plt.ylabel('Delta B')
    plt.legend()

    # Save the data
    data_to_save = {
        'time_to_nucleation': time_to_nucleation.tolist(),
        'strains': strains.tolist(),
        'Delta_Bs': Delta_Bs.tolist()
    }
    with open('time_to_nucleation.pkl', 'wb') as f:
        pickle.dump(data_to_save, f)

update_plot()
plt.draw()
plt.pause(0.001)

for i in range(len(Delta_Bs)):

    Delta_B = Delta_Bs[i]
    pfc0 = cf.PhaseFieldCrystal2DSquare(nx,nx, t=t, r=Delta_B/B_x, v=v, psi0=psi0)
    pfc0.conf_PFC_from_amplitudes()  
    pfc0.plot_lib = 'matplotlib'
    from_low_strain=True
    finished_with_Delta_B = False

    strain_counter = 0

    for m in range(len(strains)):
        if not finished_with_Delta_B:
            strain = strains[strain_counter] if from_low_strain else strains[-strain_counter-1]
            vidar_strain = 1/(1-strain) - 1
            j = strain_counter if from_low_strain else -strain_counter-1
            strain_counter += 1

            print(f"Testing Delta_B = {Delta_B:.5f} and strain = {strain:.5f}")

            pfc = copy.deepcopy(pfc0)
            distortion = np.eye(2)*vidar_strain
            pfc.conf_apply_distortion(distortion, update_q_and_a_vectors=True)
            pfc.psi = pfc.psi + 1e-2*(np.random.rand(pfc.psi.shape[0], pfc.psi.shape[1])-0.5)
            pfc.psi_f = sp.fft.fftn(pfc.psi)

            pfc.evolve_PFC(100,suppress_output=True)
            
            noise_strength = 0.01
            max_alpha = 0
            time_to_nucleation[j,i] = 0

            free_energy  = pfc.calc_free_energy()

            free_energy_sharp_drop=False
            free_energy_stable=False

            time = 0
            free_energy_old = pfc.calc_free_energy()
            initial_free_energy = free_energy_old

            trailing_free_energies = np.exp(-101/np.arange(1,102,1))*initial_free_energy

            dFdt_old = np.inf

            while not free_energy_stable and not free_energy_sharp_drop and not finished_with_Delta_B:
                pfc.evolve_PFC(100,suppress_output=True)
                time_to_nucleation[j,i] += time_step    

                free_energy_new = pfc.calc_free_energy()

                trailing_free_energies = np.roll(trailing_free_energies, -1)
                trailing_free_energies[-1] = free_energy_new

                dFdt_new = (free_energy_new - free_energy_old)/(time_step*pfc0.dt)
                d2Fdt2 = (dFdt_new - dFdt_old)/(time_step*pfc0.dt)
                dFdt_old = dFdt_new

                free_energy_old = free_energy_new

                derivatives = calc_derivatives(trailing_free_energies, time_step*pfc0.dt)

                number_of_positive_derivatives = np.sum(derivatives > 0)

                tqdm.write(f"\rTime: {time_to_nucleation[j,i]:10.0f}, dFdt: {dFdt_new:.16f}, Number of positive derivatives: {number_of_positive_derivatives:2.0f} ", end='')
                
                if number_of_positive_derivatives == 100:    
                    free_energy_stable = True
                    print("Stable: Free energy stable")
                    time_to_nucleation[i,j] = np.nan

                    if from_low_strain:
                        if strain_counter == len(strains):
                            finished_with_Delta_B = True
                            
                        from_low_strain=False
                        strain_counter=0

                    else:
                        finished_with_Delta_B = True
                    
                    update_plot()
                    plt.draw()
                    plt.pause(0.001)

                if number_of_positive_derivatives  == 0 or (initial_free_energy-free_energy_new)/abs(initial_free_energy) > 0.1:
                    free_energy_sharp_drop=True
                    print("Unstable: Free energy sharp drop")

                    update_plot()
                    plt.draw()
                    plt.pause(0.001)
                    

            free_energy = free_energy_new
        


plt.show()
# # Save the data
# with open(f"{folder_name}/time_to_nucleation.pkl", "wb") as f:
#     pickle.dump(time_to_nucleation, f)
