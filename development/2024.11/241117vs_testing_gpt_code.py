import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scipy as sp


# Initialize a Bose-Einstein Condensate (BEC) system
bec = cf.BoseEinsteinCondensate(dim=2, xRes=256, yRes=256, xmin=-20, xmax=20, ymin=-20, ymax=20, dt=0.01)

# Configure the initial condition with a vortex-antivortex pair
vortex_charge = 1  # Positive vortex
antivortex_charge = -1  # Negative vortex
vortex_position = [5, 0]  # Position of vortex
antivortex_position = [-5, 0]  # Position of antivortex
bec.conf_insert_vortex(charge=vortex_charge, position=vortex_position)
bec.conf_insert_vortex(charge=antivortex_charge, position=antivortex_position)

# Relax the system to dissipate any unwanted artifacts
bec.evolve_relax(100)

# Simulate vortex-antivortex annihilation
bec.evolve_dGPE(number_of_steps=100)

# Plot the resulting field
fig = bec.plot_complex_field(bec.psi)
fig.show()

