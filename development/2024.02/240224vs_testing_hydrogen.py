
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import matplotlib.pyplot as plt

import numpy as np
from scipy.special import sph_harm, genlaguerre
from scipy.constants import epsilon_0, hbar, e, m_e

def radial_wavefunction(n, l, r):
    """
    Calculate the radial part of the hydrogen atom wavefunction.
    
    Parameters:
    - n: Principal quantum number
    - l: Orbital quantum number
    - r: Radial distance
    
    Returns:
    - R_nl: The radial part of the wavefunction
    """
    # Normalization constant
    a0 = 4 * np.pi * epsilon_0 * hbar**2 / (m_e * e**2)  # Bohr radius
    rho = 2 * r / (n * a0)
    normalization = np.sqrt((2 / (n * a0))**3 * np.math.factorial(n - l - 1) / (2 * n * np.math.factorial(n + l)))
    # Radial part
    L_nl = genlaguerre(n-l-1, 2*l+1)  # Associated Laguerre polynomial
    R_nl = normalization * np.exp(-rho / 2) * rho**l * L_nl(rho)
    return R_nl

def angular_wavefunction(l, m, theta, phi):
    """
    Calculate the angular part of the hydrogen atom wavefunction.
    
    Parameters:
    - l: Orbital quantum number
    - m: Magnetic quantum number
    - theta: Colatitude (angle from the z-axis)
    - phi: Azimuthal angle (angle from the x-axis in the xy-plane)
    
    Returns:
    - Y_lm: The angular part of the wavefunction (spherical harmonics)
    """
    Y_lm = sph_harm(m, l, phi, theta)
    return Y_lm

def hydrogen_wavefunction(n, l, m, r, theta, phi):
    """
    Calculate the hydrogen atom wavefunction.
    
    Parameters:
    - n: Principal quantum number
    - l: Orbital quantum number
    - m: Magnetic quantum number
    - r: Radial distance
    - theta: Colatitude
    - phi: Azimuthal angle
    
    Returns:
    - psi: The hydrogen atom wavefunction
    """
    R_nl = radial_wavefunction(n, l, r)
    Y_lm = angular_wavefunction(l, m, theta, phi)
    psi = R_nl * Y_lm
    return psi


qm = cf.QuantumMechanics(3,xlim=[-2,2],ylim=[-2,2],zlim=[-2,2])

r = np.sqrt(qm.x**2 + qm.y**2 + qm.z**2)
theta = np.arccos(qm.z / r)
phi = np.arctan2(qm.y, qm.x)

qm.psi = hydrogen_wavefunction(1, 0, 0, r, theta, phi)
print(qm.psi)
qm.plot_complex_field(qm.psi)

plt.show()