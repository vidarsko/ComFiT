import pickle
import matplotlib.pyplot as plt
import numpy as np

# Path to the pickle file
file_path = 'time_to_nucleation.pkl'

# Open and read the pickle file
with open(file_path, 'rb') as file:
    data = pickle.load(file)

    time_to_nucleation = np.array(data['time_to_nucleation'])
    strains = data['strains']
    Delta_Bs = data['Delta_Bs']

Delta_Bs_mesh, strains_mesh = np.meshgrid(Delta_Bs,strains)

data_file_path = 'C:/Users/vidar/Desktop/ComFiT/development/2024.10/stab7mode_98.dat'
data = np.loadtxt(data_file_path)
# Save column one into strain_ana and column two into deltab_ana
strain_ana = data[:, 0]
deltab_ana = data[:, 1]


plt.pcolormesh(strains_mesh.transpose(), Delta_Bs_mesh.transpose(), np.log10(time_to_nucleation.transpose()), shading='auto')
plt.colorbar(label='Time to Nucleation (powers of 10)')
plt.plot(strain_ana, deltab_ana, 'r', label='Analytical')
plt.xlabel('Strain')
plt.ylabel('Delta B')
plt.legend()

plt.show()