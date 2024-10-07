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

nx=40
B_x=0.98
psi0=0
r = 0.02/B_x #Delta B
t=1/2/B_x
v=1/3/B_x
time_limit = 100000
time_step = 1000

# Create a folder with the current date and time
folder_name = datetime.now().strftime("%Y%m%d_%H%M%S")+f"_nx_{nx}_Bx_{B_x}_psi0_{psi0}_r_{B_x*r}_t_{B_x*t}_v_{B_x*v}_time_limit_{time_limit}"
if not os.path.exists(folder_name):
    os.makedirs(folder_name)


# pfc.plot_lib = 'matplotlib'

res = 101
Delta_Bs = np.linspace(-0.02,0.06,res)
strains = np.linspace(-0.1,0.1,res)

# delta_Blimit = (B_x*t)**2/(15*B_x*v) - B_x*(1-(1-strains)**2)**2
# plt.plot(strains,delta_Blimit)
# plt.show()
# raise Exception("This script is not finished yet. Please do not run it.")
Delta_Bs_mesh, strains_mesh = np.meshgrid(Delta_Bs,strains)
time_to_nucleation = np.nan*np.zeros_like(Delta_Bs_mesh)


def plot():
    fig, axs = plt.subplots(1,2)
    pfc.plot_PFC(ax=axs[0],title='PFC')

    alpha = pfc.calc_dislocation_density()
    pfc.plot_field(np.sqrt(alpha[0]**2+alpha[1]**2), vlim=[0,0.0045], ax=axs[1], title='Dislocation density')

for i in range(len(Delta_Bs)):
    Delta_B = Delta_Bs[i]
    pfc0 = cf.PhaseFieldCrystal2DTriangular(nx,round(nx/np.sqrt(3)), t=t, r=Delta_B/B_x, v=v, psi0=psi0)
    pfc0.conf_PFC_from_amplitudes()  
    pfc0.plot_lib = 'matplotlib'
    from_low_strain=True
    finished_with_Delta_B = False
    counter = 0
    for m in range(len(strains)):
        if not finished_with_Delta_B:
            strain = strains[counter] if from_low_strain else strains[-counter-1]
            j = counter if from_low_strain else -counter-1
            counter += 1
            print("Testing Delta_B = ", Delta_B, " and strain = ", strain)

            pfc = copy.deepcopy(pfc0)
            distortion = np.eye(2)*strain
            pfc.conf_apply_distortion(distortion)
            
            noise_strength = 0.01
            max_alpha = 0
            time_to_nucleation[i,j] = 0
            
            time_array = np.arange(0,time_limit,time_step)

            A1max = np.nan*np.zeros(time_array.shape)
            A1mean = np.nan*np.zeros(time_array.shape)

            A2max = np.nan*np.zeros(time_array.shape)
            A2mean = np.nan*np.zeros(time_array.shape)

            A3max = np.nan*np.zeros(time_array.shape)
            A3mean = np.nan*np.zeros(time_array.shape)

            is_liquid = False

            counter=0

            def update_variables():
                alpha = pfc.calc_dislocation_density()
                max_alpha = np.sqrt(alpha[0]**2+alpha[1]**2).max()
                
                eta = pfc.calc_demodulate_PFC()
                eta0 = np.abs(eta[0])
                A1max[counter] = eta0.max()
                A1mean[counter] = eta0.mean()

                eta1 = np.abs(eta[1])
                A2max[counter] = eta1.max()
                A2mean[counter] = eta1.mean()

                eta2 = np.abs(eta[2])
                A3max[counter] = eta2.max()
                A3mean[counter] = eta2.mean()

                max_amp = np.max([A1max[int(time_to_nucleation[i,j]/time_step)], A2max[int(time_to_nucleation[i,j]/time_step)], A3max[int(time_to_nucleation[i,j]/time_step)]])
                return max_alpha, max_amp

            max_alpha, max_amp = update_variables()

            if Delta_B > (B_x*pfc.t)**2/(15*B_x*pfc.v) - B_x*(1-(1-strain)**2)**2:
                print("Delta_B too high compared to strain - is liquid")
                is_liquid = True
                time_to_nucleation[i,j] = np.nan

            while max_alpha < 0.004 and time_to_nucleation[i,j] < time_limit and not is_liquid:
                pfc.psi = pfc.psi + noise_strength*np.random.randn(*pfc.psi.shape)
                pfc.psi_f = sp.fft.fftn(pfc.psi)
                pfc.evolve_PFC(time_step)
                counter += 1

                time_to_nucleation[i,j] += time_step

                max_alpha, max_amp = update_variables()

                print(max_amp)
                if max_amp < 0.001:
                    is_liquid = True
                    time_to_nucleation[i,j] = np.nan
                # print("Max alpha = ", max_alpha)
                # plt.pause(0.1)
            
            if time_to_nucleation[i,j] == time_limit:
                if from_low_strain:
                    from_low_strain = False
                    counter=0
                else:
                    finished_with_Delta_B = True

            # plt.clf()
            # plot()
            time = 'nan' if np.isnan(time_to_nucleation[i,j]) else int(time_to_nucleation[i,j])
            with open(os.path.join(folder_name, f"pfc_DeltaB_{Delta_B:+07.5f}_strain_{strain:+01.5f}_time_{time}.pkl"), 'wb') as pkl_file:
                pickle.dump(pfc, pkl_file)

            data = {
                'Delta_B': Delta_B,
                'strain': strain,
                'time_to_nucleation': time_to_nucleation[i,j],
                'A1max': A1max.tolist(),
                'A1mean': A1mean.tolist(),
                'A2max': A2max.tolist(),
                'A2mean': A2mean.tolist(),
                'A3max': A3max.tolist(),
                'A3mean': A3mean.tolist()
            }

            filename = os.path.join(folder_name, f"data_DeltaB_{Delta_B:+07.5f}_strain_{strain:.2f}_time_{time}.json")
            with open(filename, 'w') as f:
                json.dump(data, f)
            # plt.savefig(os.path.join(folder_name, f"DeltaB_{Delta_B:+07.5f}_strain_{strain:.2f}_time_{time}.png"))
            # pfc.plot_PFC()
            # plt.show()

            # plt.clf()
            # plt.pcolormesh(strains_mesh, Delta_Bs_mesh, time_to_nucleation.transpose(), shading='auto')
            # plt.colorbar(label='Time to Nucleation')
            # plt.savefig(os.path.join(folder_name, "time_to_nucleation.png"))
            # plt.show()
        



