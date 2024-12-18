import comfit as cf
import os
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import copy
import pickle
from tqdm import tqdm

import cupy as cp   

def create_fft_functions(pfc):

    def fft(field):
        if not isinstance(field, cp.ndarray):
            field = cp.asarray(field)
        return cp.asnumpy(cp.fft.fftn(field, axes=range(-pfc.dim, 0)))

    def ifft(field):
        if not isinstance(field, cp.ndarray):
            field = cp.asarray(field)
        return cp.asnumpy(cp.fft.ifftn(field, axes=range(-pfc.dim, 0)))

    return fft, ifft

nx=20
B_x=0.98
psi0=0
t=1/2/B_x
v=1/3/B_x

nsteps=100

Delta_Bs = np.arange(-0.15, 0.07, 0.005)
strains = np.arange(-0.2, 0.205, 0.005)

Delta_Bs_mesh, strains_mesh = np.meshgrid(Delta_Bs,strains)
time_to_nucleation = np.nan*np.zeros_like(Delta_Bs_mesh)

strain_type = 'shear'
pfc_type = 'square'

current_time = datetime.now().strftime(f"%y%m%d - Simulations %H%M%S - {pfc_type} - {strain_type}")
results_folder = os.path.join(r'C:\Users\vidar\UiO Dropbox\Vidar Skogvoll\My projects', current_time)
if not os.path.exists(results_folder):
    os.makedirs(results_folder)

def update_plot():
    plt.clf()
    plt.pcolormesh(strains_mesh.transpose(), Delta_Bs_mesh.transpose(), np.log10(time_to_nucleation.transpose()), shading='auto')
    plt.colorbar(label='Time to Nucleation (powers of 10)')
    plt.xlabel('Strain')
    plt.ylabel('Delta B')
    plt.legend()

    plt.savefig(os.path.join(results_folder, 'nucleation_plot.png'))

    # Save the data
    data_to_save = {
        'time_to_nucleation': time_to_nucleation.tolist(),
        'strains': strains.tolist(),
        'Delta_Bs': Delta_Bs.tolist()
    }
    with open(os.path.join(results_folder, 'time_to_nucleation.pkl'), 'wb') as f:
        pickle.dump(data_to_save, f)

    # plt.draw()
    # plt.pause(0.001)

update_plot()


def save_PFC_plot(Delta_B, strain, pfc):
    plt.clf()
    pfc.plot_field(pfc.psi)
    plt.title(f'Delta B = {Delta_B:.5f}, Strain = {strain:.5f}')
    from_type_strain = 'low' if from_low_strain else 'high'
    plt.savefig(os.path.join(results_folder, f'PFC_Delta_B_{Delta_B:.5f}_{from_type_strain}.png'))


for i in range(len(Delta_Bs)):

    Delta_B = Delta_Bs[i]
    pfc0 = cf.PhaseFieldCrystal2DSquare(nx,nx, t=t, r=Delta_B/B_x, v=v, psi0=psi0)
    pfc0.conf_PFC_from_amplitudes()  
    pfc0.plot_lib = 'matplotlib'

    # Enable gpu acceleration
    # pfc0.fft, pfc0.ifft = create_fft_functions(pfc0)

    from_low_strain=True
    finished_with_Delta_B = False

    strain_counter = 0

    for m in range(len(strains)):
        if not finished_with_Delta_B:
            strain = strains[strain_counter] if from_low_strain else strains[-strain_counter-1]
            # vidar_strain = 1/(1-strain) - 1
            vidar_strain = strain
            j = strain_counter if from_low_strain else -strain_counter-1
            strain_counter += 1

            print(f"Testing Delta_B = {Delta_B:.5f}, strain = {strain:.5f}, pfc_type = {pfc_type}, strain_type = {strain_type}")

            pfc = copy.deepcopy(pfc0)
            distortion = np.array([[0,vidar_strain],[vidar_strain,0]])
            pfc.conf_apply_distortion(distortion, update_q_and_a_vectors=True)
            pfc.psi = pfc.psi + 1e-2*(np.random.rand(pfc.psi.shape[0], pfc.psi.shape[1])-0.5)
            pfc.psi_f = sp.fft.fftn(pfc.psi)

            pfc.evolve_PFC(50,suppress_output=True)
            
            time_to_nucleation[j,i] = 0
            free_energy  = pfc.calc_free_energy()

            free_energy_sharp_drop=False
            free_energy_stable=False

            free_energy_old = pfc.calc_free_energy()
            initial_free_energy = free_energy_old

            trailing_free_energies = np.exp(-101/np.arange(1,102,1))*initial_free_energy

            dFdt_old = np.inf

            while not free_energy_stable and not free_energy_sharp_drop and not finished_with_Delta_B:
                pfc.evolve_PFC(nsteps,suppress_output=True)
                time_to_nucleation[j,i] += nsteps    

                free_energy_new = pfc.calc_free_energy()
                dFdt = (free_energy_new - free_energy_old)/(nsteps*pfc.dt)
                free_energy_old = free_energy_new

                tqdm.write(f"\rTime: {time_to_nucleation[j,i]:10.0f}, dFdt: {dFdt:.16f} ", end='')
                
                if abs(dFdt) < 1e-10:    
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


                if (initial_free_energy-free_energy_new)/abs(initial_free_energy) > 0.1:
                    free_energy_sharp_drop=True
                    print("Unstable: Free energy sharp drop")

                    save_PFC_plot(Delta_B, strain, pfc)

                    update_plot()
                    
                    

            free_energy = free_energy_new
        

plt.show()