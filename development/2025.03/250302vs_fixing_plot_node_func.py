import comfit as cf
import matplotlib.pyplot as plt
import numpy as np

### Task 1: initialise a bec in two dimensions with resolution 101 in x and y. Make gamma = 0.05 
bec = cf.BoseEinsteinCondensate(2,xRes=101,yRes=101,gamma=0.05,dt=0.1, plot_lib = 'matplotlib')

### First we set the size of the harmonic
R_tf = 40



### Here we set the size and velocity of the stirrer
stirrer_radius = 20
stirrer_velocity = 0.6
freq = stirrer_velocity/stirrer_radius
size =4
strength = .9

### Defining the function for the time-dependent potential
def V_t(t):
    pos_x = bec.xmid + stirrer_radius * np.cos(freq * t)
    pos_y = bec.ymid + stirrer_radius * np.sin(freq * t)
    stirrer = bec.calc_Gaussian(width=size/np.sqrt(2), top=strength, position=[pos_x, pos_y])
    harmonic = bec.calc_harmonic_potential(R_tf)
    return   harmonic + stirrer



### Task 2: Set the potential to the above function, initialise the Thomas Fermi 
# ground state and relax the system using the  evolve_relax(...) solver for 20 time steps

bec.conf_external_potential(V_t, additive=False)

bec.conf_initial_condition_Thomas_Fermi()

bec.evolve_relax(50) 

# bec.plot_field(np.abs(bec.psi)**2,cmap_symmetric=False,colormap = 'viridis')

### Task 4: Evolve the system with the time-dependent potential using the ETD4RK scheme

bec.evolve_dGPE( 300, method='ETD4RK') 


### Task 5: Track the defects and their velocity and plot the result 

nodes = bec.calc_vortex_nodes()

fig,ax=bec.plot_field(np.abs(bec.psi)**2,cmap_symmetric=False, colormap = 'gray')
fig, ax = bec.plot_nodes(nodes, fig=fig, ax=ax)
bec.show(fig)