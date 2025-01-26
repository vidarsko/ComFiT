import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scipy as sp
import pickle
import glob


res = 101
Delta_Bs = np.linspace(-0.02,0.06,res)
strains = np.linspace(-0.1,0.1,res)

delta_B = Delta_Bs[70]
strain = strains[70]

print(f"Strain: {strain}, Delta_B: {delta_B}")

folder_name = r'C:\Users\vidar\UiO Dropbox\Vidar Skogvoll\My projects\241001 - Data from exploring\nx_20_Bx_0.98_psi0_0_r_0.02_t_0.5_v_0.3333333333333333_time_limit_100000_20241002_230447'

# Find the file that matches the pattern
pattern = os.path.join(folder_name, f"pfc_DeltaB_{delta_B:+07.5f}_strain_{strain:+01.5f}_time_*.pkl")
files = glob.glob(pattern)

if len(files) != 1:
    raise FileNotFoundError(f"Expected one file, but found {len(files)} files.")

file_path = files[0]

# Load the data from the file
with open(file_path, 'rb') as file:
    data = pickle.load(file)

pfc = data
# pfc.plot_lib = 'matplotlib'
pfc.plot_field(pfc.psi)
plt.show()