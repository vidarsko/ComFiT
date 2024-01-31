import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

# Testing the wave equation in two dimensions
# Initialize BaseSystem object
bs = cf.BaseSystem(2, xRes=51, dx=1, yRes=51, dy=1)

# Add linear and non-linear functions
def calc_omega_f(bs):
    return -bs.calc_k2()

def calc_nonlinear_evolution_function_f(bs,T):
    f = bs.A*(bs.T0-T)*np.exp(-((bs.x-bs.xmid)**2+(bs.y-bs.ymid)**2)/(2*bs.sigma**2))
    return sp.fft.fft2(f)

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
    T = T_flat.reshape((bs.xRes, bs.yRes))
    forcing = bs.A * (bs.T0 - T) * np.exp(-((bs.x - bs.xmid) ** 2 + (bs.y - bs.ymid) ** 2) / (2 * bs.sigma ** 2))
    return forcing.flatten()

# Function representing the discretized PDE
def heat_equation(t,T_flat):
    T = T_flat.reshape((bs.xRes, bs.yRes))
    dTdx2 = np.roll(T, -1, axis=0) - 2 * T + np.roll(T, 1, axis=0)
    dTdx2 /= bs.dx ** 2
    dTdy2 = np.roll(T, -1, axis=1) - 2 * T + np.roll(T, 1, axis=1)
    dTdy2 /= bs.dy ** 2
    # print("dTdx2.flatten shape:", dTdx2.flatten().shape)
    # print("dTdy2.flatten shape:", dTdy2.flatten().shape)
    # print("forcing_term(T_flat) shape:", forcing_term(T_flat).shape)
    return dTdx2.flatten() + dTdy2.flatten() + forcing_term(T_flat)


# Set test parameters
bs.A = 0.1
bs.sigma = 5
bs.T0 = 20

# ComFiT simulation

# Initial condition
bs.T = np.zeros((bs.xRes,bs.yRes))
bs.T_f = sp.fft.fft2(bs.T)

# Scipy benchmark

# Initial condition
T_initial = np.zeros((bs.xRes,bs.yRes))

fig, axs = plt.subplots(1,2,figsize=(10,5))

colorbar_on = True
for n in range(51):
    print(n)
    # Evolve the system
    bs.evolve(100)
    

    axs[0].cla()
    bs.plot_field(bs.T,ax=axs[0],colormap='bluewhitered',cmap_symmetric=False,
                clims=[0,bs.T0], colorbar=colorbar_on)
    axs[0].set_title(f'ComFiT, method=ETD2RK,\n dt={bs.dt},dx={bs.dx},dy = {bs.dy}')

    # Time span for the integration
    t_span = (0, 10)

    # Solve the equation
    T_benchmark = sp.integrate.solve_ivp(heat_equation, t_span, T_initial.flatten(), method='RK45')
    T_initial = T_benchmark.y[:, -1].reshape((bs.xRes, bs.yRes))

    axs[1].cla()
    bs.plot_field(T_initial,ax=axs[1],colormap='bluewhitered',cmap_symmetric=False,
                clims=[0,bs.T0],colorbar=colorbar_on)
    axs[1].set_title(f'Scipy.integrate.solve_ivp,\n method=RK45, dx={bs.dx}, dy = {bs.dy}')
    
    fig.suptitle(f't={n*10}, sigma={bs.sigma}, A={bs.A}, T0={bs.T0}')

    cf.tool_save_plot(n)

    colorbar_on = False
cf.tool_make_animation_gif(n)

