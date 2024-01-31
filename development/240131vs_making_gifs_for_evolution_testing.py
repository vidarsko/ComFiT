import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

# Testing the wave equation in one dimensions
# Initialize BaseSystem object
bs = cf.BaseSystem(1, xRes=101, dx=1)

# Add linear and non-linear functions
def calc_omega_f(bs):
    return -bs.calc_k2()

def calc_nonlinear_evolution_function_f(bs,T):
    f = bs.A*(bs.T0-T)*np.exp(-(bs.x-bs.xmid)**2/(2*bs.sigma**2))
    return sp.fft.fft(f)

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
def forcing_term(T):
    return bs.A*(bs.T0-T)*np.exp(-(bs.x-bs.xmid)**2/(2*bs.sigma**2))

# Function representing the discretized PDE
def heat_equation(t, T):
    dTdx2 = np.roll(T, -1) - 2 * T + np.roll(T, 1)
    dTdx2 /= bs.dx**2
    return dTdx2 + forcing_term(T)


# Set test parameters
bs.A = 0.1
bs.sigma = 5
bs.T0 = 20

# ComFiT simulation

# Initial condition
bs.T = np.zeros((bs.xRes))
bs.T_f = sp.fft.fft(bs.T)   

# Scipy benchmark

# Initial condition
T_initial = np.zeros((bs.xRes))

fig, axs = plt.subplots(1,2,figsize=(10,5))
for n in range(10):
    # Evolve the system
    bs.evolve(10)

    axs[0].cla()
    axs[0].plot(bs.x,bs.T)
    axs[0].set_ylim([0,bs.T0])
    axs[0].set_title('ComFiT')
    axs[0].set_xlabel('x')
    axs[0].set_ylabel('T')
    axs[0].grid()

    # Time span for the integration
    t_span = (0, 1)

    # Solve the equation
    T_benchmark = sp.integrate.solve_ivp(heat_equation, t_span, T_initial, method='RK45')

    T_initial = T_benchmark.y[:,-1]
    axs[1].cla()
    axs[1].plot(bs.x,T_initial)
    axs[1].set_ylim([0,bs.T0])
    axs[1].set_title('Scipy')
    axs[1].set_xlabel('x')
    # axs[1].set_ylabel('T')
    axs[1].grid()
    
    cf.tool_save_plot(n)
cf.tool_make_animation_gif(n)

