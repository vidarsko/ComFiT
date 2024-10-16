import os
import numpy as np
import matplotlib.pyplot as plt
import json
from datetime import datetime

# Import data from 241003vs_strainfit.dats file
data_file_path = 'C:/Users/vidar/Desktop/ComFiT/development/2024.10/stab7mode_98.dat'
# data_file_path = 'C:/Users/vidar/Desktop/ComFiT/development/2024.10/241003vs_strainfit.dats'
# data_file_path = 'C:/Users/vidar/UiO Dropbox/Vidar Skogvoll/My projects/241001 - Data from exploring/nx_20_Bx_0.98_psi0_0_r_0.02_t_0.5_v_0.3333333333333333_time_limit_100000_20241002_230447'
data = np.loadtxt(data_file_path)

# Save column one into strain_ana and column two into deltab_ana
strain_ana = data[:, 0]
deltab_ana = data[:, 1]

# Parameters
nx = 20
B_x = 0.98
psi0 = 0
r = 0.02 / B_x  # Delta B
t = 1 / 2 / B_x
v = 1 / 3 / B_x
time_limit = 100000

# Create a folder with the current date and time
# folder_name = 'nx_20_Bx_0.98_psi0_0_r_0.02_t_0.5_v_0.3333333333333333_time_limit_100000_20241002_230447'
folder_name = r'C:\Users\vidar\Desktop\ComFiT\20241016_213811_nx_40_Bx_0.98_psi0_0_t_0.5_v_0.3333333333333333_time_limit_100000'
# if not os.path.exists(folder_name):
#     os.makedirs(folder_name)

# Parameters for plotting
res = 101
Delta_Bs = np.linspace(-0.15,0.06,res)
strains = np.linspace(-0.2,0.2,res)
Delta_Bs_mesh, strains_mesh = np.meshgrid(Delta_Bs, strains)
time_to_nucleation = np.nan * np.zeros_like(Delta_Bs_mesh)

# Read JSON files and extract data
for json_file in os.listdir(folder_name):
    if json_file.endswith(".json"):
        with open(os.path.join(folder_name, json_file), 'r') as f:
            data = json.load(f)
            Delta_B = data['Delta_B']
            strain = data['strain']
            time_to_nucleation_value = data['time_to_nucleation']
            
            # Find the closest indices in the mesh grid
            i = (np.abs(Delta_Bs - Delta_B)).argmin()
            j = (np.abs(strains - strain)).argmin()
            
            time_to_nucleation[i, j] = time_to_nucleation_value

# Plotting
plt.figure()
plt.pcolormesh(strains_mesh, Delta_Bs_mesh, time_to_nucleation.transpose(), shading='auto')
plt.colorbar(label='Time to Nucleation')
plt.xlabel('Strain')
plt.ylabel('Delta B')
plt.title('Time to Nucleation as a function of Strain and Delta B')
plt.plot(strain_ana, deltab_ana, 'r', label='Analytical')
plt.savefig(os.path.join(folder_name, "time_to_nucleation.png"))
plt.show()
