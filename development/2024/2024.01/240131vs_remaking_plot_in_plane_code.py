import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

# Testing the wave equation in two dimensions
# Initialize BaseSystem object
bs = cf.BaseSystem(3, xRes=21, dx=1, yRes=21, dy=1, zRes=21, dz=1)

# Add linear and non-linear functions
def calc_omega_f(bs):
    return -bs.calc_k2()

def calc_nonlinear_evolution_function_f(bs,T):
    f = bs.A*(bs.T0-T)*np.exp(-((bs.x-bs.xmid)**2+(bs.y-bs.ymid)**2+(bs.z-bs.zmid)**2)/(2*bs.sigma**2))
    return sp.fft.fftn(f)

# Add functions to BaseSystem object
bs.calc_omega_f = calc_omega_f.__get__(bs)
bs.calc_nonlinear_evolution_function_f = calc_nonlinear_evolution_function_f.__get__(bs)

# Define evolve function
def evolve(bs,number_of_steps, method = 'ETD2RK'):
    omega_f = bs.calc_omega_f()

    integrating_factors_f, solver = bs.calc_integrating_factors_f_and_solver(omega_f,method)

    for n in range(number_of_steps):
        bs.T, bs.T_f = solver(integrating_factors_f,
                                bs.calc_nonlinear_evolution_function_f,
                                bs.T, bs.T_f)
        bs.T = np.real(bs.T)

# Add evolve function to BaseSystem object
bs.evolve = evolve.__get__(bs)

# Scipy benchmark functions

# Forcing term function
def forcing_term(T_flat):
    T = T_flat.reshape((bs.xRes, bs.yRes, bs.zRes))
    forcing = bs.A * (bs.T0 - T) * np.exp(-((bs.x - bs.xmid) ** 2 + (bs.y - bs.ymid) ** 2 + (bs.z - bs.zmid) ** 2) / (2 * bs.sigma ** 2))
    return forcing.flatten()

# Function representing the discretized PDE
def heat_equation(t,T_flat):
    T = T_flat.reshape((bs.xRes, bs.yRes, bs.zRes))
    dTdx2 = np.roll(T, -1, axis=0) - 2 * T + np.roll(T, 1, axis=0)
    dTdx2 /= bs.dx ** 2
    dTdy2 = np.roll(T, -1, axis=1) - 2 * T + np.roll(T, 1, axis=1)
    dTdy2 /= bs.dy ** 2
    dTdz2 = np.roll(T, -1, axis=2) - 2 * T + np.roll(T, 1, axis=2)
    dTdz2 /= bs.dz ** 2
    return dTdx2.flatten() + dTdy2.flatten() + dTdz2.flatten() + forcing_term(T_flat)

# Set test parameters
bs.A = 0.1
bs.sigma = 3
bs.T0 = 20

# ComFiT simulation

# Initial condition
bs.T = np.zeros((bs.xRes,bs.yRes,bs.zRes))
bs.T_f = sp.fft.fftn(bs.T)

# Scipy benchmark

# Initial condition
T_initial = np.zeros((bs.xRes,bs.yRes,bs.zRes))


colorbar_on = True


# Evolve the system
bs.evolve(1000)


# ax = bs.plot_field_in_plane2(bs.T,normal_vector=[1,-1,0], colormap='bluewhitered',clims=[0,bs.T0])
# ax = bs.plot_field(bs.T,ax=ax,colormap='bluewhitered',cmap_symmetric=False,
            # clims=[0,bs.T0], colorbar=colorbar_on, number_of_layers=5)
# ax.set_title(f'ComFiT, method=ETD2RK,\n dt={bs.dt},dx={bs.dx},dy = {bs.dy}')

bs.plot_field_in_plane_with_mayavi(bs.T,normal_vector=[1,-1,0], colormap='bluewhitered')

# cf.tool_save_plot_matplotlib(n)

plt.show()
# cf.tool_make_animation_gif(n)

