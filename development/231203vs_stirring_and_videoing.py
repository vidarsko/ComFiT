import comfit as cf
import matplotlib.pyplot as plt
import numpy as np

bec = cf.BEC(2,xRes=100,yRes=100,gamma=0.05,dt=0.1)

### First we set the size of the harmonic
R_tf = 40

### Here we set the size and velocity of the stirrer
stirrer_radius = 20
stirrer_velocity = 0.6
freq = stirrer_velocity/stirrer_radius
size = 4
strength = .9

### Defining the function for the time-dependent potential
def V_t():
    pos_x = bec.xmid + stirrer_radius * np.cos(freq * bec.t)
    pos_y = bec.ymid + stirrer_radius * np.sin(freq * bec.t)
    stirrer = bec.calc_gaussian_stirring_potential(size, strength, [pos_x, pos_y])
    harmonic = bec.set_harmonic_potential(R_tf)
    return   harmonic + stirrer

### Set the potential to the t=0 value, initialise the TF ground-state and relax in imaginary time

bec.V0 = V_t()
bec.conf_initial_condition_Thomas_Fermi()
bec.evolve_relax_BEC(20,'ETD2RK')

### Updating the potential to the time-dependent function V_t()
bec.set_time_dependent_potential(V_t)




for n in range(1000):
    plt.gcf().clf()
    bec.evolve_dGPE( 20,'ETD2RK')
    psi_old = bec.psi
    bec.evolve_dGPE(1)
    bec.plot_field(np.abs(bec.psi)**2,cmap_symmetric=False,colormap = 'winter',colorbar=False)
    Dnodes = bec.calc_vortex_nodes((bec.psi-psi_old)/bec.dt)
    bec.plot_vortex_nodes(Dnodes)
    cf.tool_save_plot(n)
cf.tool_make_animation(n)