import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scipy as sp


bec = cf.BoseEinsteinCondensate(2,xRes=256,yRes=128,gamma=0,dt=0.05)

### task 1: Set the potiential to a constant gaussian at this position [bec.xmid+50,bec.ymid] with size = 5 and
# strength = 4
#Set the potential by defininf V0 = ...
# and calling the function bec.conf_external_potential(V0)

bec = cf.BoseEinsteinCondensate(2,xRes=256,yRes=128,gamma=0,dt=0.05)

pot = bec.calc_gaussian_stirring_potential(5, 4, [bec.xmid+50,bec.ymid] )
bec.conf_external_potential(pot, additive=False)

### task 2: Initialise the wavefunction using the Thomas-Fermi ground state and relax the system in imaginary time 
# for 50 time steps. Plot the absolut value squared of the wave function
#(Hint: the evolvers vere discussed in the previous notebook)

bec.conf_initial_condition_Thomas_Fermi()

bec.evolve_relax(50, method='ETD2RK') 

bec.conf_dissipative_frame(interface_width=7, frame_width_x=200, frame_width_y=100)

### task 4. evolve the system in the comoving frame with vel_x = 0.4. Make an animation of the absolut value squared 
# of the wavefunction and mark the position of the defects.   



vel_x = 0.40



bec.evolve_comoving_dGPE(1000,vel_x,method='ETD4RK')

psi_prev = np.copy(bec.psi)

bec.evolve_comoving_dGPE(10,vel_x,method='ETD4RK')

dt_psi = (bec.psi - psi_prev)/(bec.dt*10)
nodes = bec.calc_vortex_nodes(dt_psi=dt_psi)

fig=bec.plot_field(np.abs(bec.psi)**2,vlim=[0,1],colormap = 'gray')
bec.plot_vortex_nodes(nodes, fig=fig)
fig.show()
