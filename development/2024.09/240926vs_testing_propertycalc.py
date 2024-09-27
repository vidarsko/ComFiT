import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scipy as sp
from comfit.tool.tool_print_in_color import tool_print_in_color

from comfit.plot.plot_field_matplotlib import plot_field_matplotlib

# pfc = cf.PhaseFieldCrystal2DTriangular(20,round(20/np.sqrt(3)))
pfc = cf.PhaseFieldCrystal2DTriangular(1,1)
pfc.conf_PFC_from_amplitudes()
# pfc = cf.PhaseFieldCrystal2DSquare(20,20)
print(pfc.psi0)
print(pfc.A)

# Evolve PFC to reach initial state of equilibrium (100 time steps ok)
number_of_initial_steps = 100
pfc.evolve_PFC(number_of_initial_steps, suppress_output=True)
number_of_steps = 10

avg_free_energy = pfc.calc_free_energy()/(pfc.xmax*pfc.ymax)

strain = 0
# Calculate free energy to compare
def print_and_show_state(free_energy):
    # free_energy_density,_ = pfc.calc_PFC_free_energy_density_and_chemical_potential()
    # print('Sum of free energy densities: ', np.sum(free_energy_density))
    # print('Average free energy at strain: ', strain, ' is: ', avg_free_energy)
    # print('Value of dx:', pfc.dx)
    # plot_field_matplotlib(pfc, pfc.psi)
    # plt.show()
    pass

print_and_show_state(avg_free_energy)


strain_increment = 0.000001
strain += strain_increment

# Unaltered k-vectors
k0 = pfc.k[0].copy()
dx0 = pfc.dx
x0 = pfc.x.copy()
xmax0 = pfc.xmax

if pfc.dim > 1:
    k1 = pfc.k[1].copy()
    dy0 = pfc.dy
    y0 = pfc.y.copy()
    ymax0 = pfc.ymax
    
if pfc.dim > 2:
    k2 = pfc.k[2].copy()
    dz0 = pfc.dz
    z0 = pfc.z.copy()
    zmax0 = pfc.zmax

def update_lengths(pfc, strain):
    pfc.k[0] = k0/(1+strain)
    pfc.dx = dx0*(1+strain)
    pfc.x = x0*(1+strain)
    pfc.xmax = xmax0*(1+strain)
    
    if pfc.dim > 1:
        pfc.k[1] = k1/(1+strain)
        pfc.dy = dy0*(1+strain)
        pfc.y = y0*(1+strain)
        pfc.ymax = ymax0*(1+strain)
        
    if pfc.dim > 2:
        pfc.k[2] = k2/(1+strain)
        pfc.dz = dz0*(1+strain)
        pfc.z = z0*(1+strain)
        pfc.zmax = zmax0*(1+strain)

    pfc.dV = pfc.dx*pfc.dy*pfc.dz

update_lengths(pfc,strain)

# Evolve PFC
pfc.evolve_PFC(number_of_steps)
avg_free_energy_tmp = pfc.calc_free_energy()/(pfc.xmax*pfc.ymax)
print_and_show_state(avg_free_energy_tmp)

positive_strain_required = False
while avg_free_energy_tmp < avg_free_energy:
    positive_strain_required = True

    avg_free_energy = avg_free_energy_tmp

    strain += strain_increment

    update_lengths(pfc,strain)
    
    pfc.evolve_PFC(number_of_steps, suppress_output=True)
    avg_free_energy_tmp = pfc.calc_free_energy()/(pfc.xmax*pfc.ymax)
    print_and_show_state(avg_free_energy_tmp)

if positive_strain_required:
    # Going one back to get the lowest free energy
    final_strain = strain - strain_increment
    update_lengths(pfc,strain)
    tool_print_in_color('Lowest average free energy found at strain: ' + str(final_strain), 'green')

else: #negative strain required

    strain = - strain_increment

    update_lengths(pfc,strain)
    pfc.evolve_PFC(number_of_steps, suppress_output=True)
    avg_free_energy_tmp = pfc.calc_free_energy()/(pfc.xmax*pfc.ymax)
    print_and_show_state(avg_free_energy_tmp)

    while avg_free_energy_tmp < avg_free_energy:

        avg_free_energy = avg_free_energy_tmp

        strain -= strain_increment

        update_lengths(pfc, strain)
        avg_free_energy_tmp = pfc.calc_free_energy()/(pfc.xmax*pfc.ymax)
        print_and_show_state(avg_free_energy_tmp)

    # Going one back to get the lowest free energy
    final_strain = strain + strain_increment
    update_lengths(pfc, final_strain)
    tool_print_in_color('Lowest average free energy found at strain: ' + str(final_strain), 'green')   