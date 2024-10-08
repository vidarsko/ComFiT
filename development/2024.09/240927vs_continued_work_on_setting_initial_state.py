import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scipy as sp
from comfit.tool.tool_print_in_color import tool_print_in_color
import copy

from scipy.optimize import curve_fit

# pfc = cf.PhaseFieldCrystal2DTriangular(1,1, for_properties_calculation=True)
pfc = cf.PhaseFieldCrystal2DSquare(1,1, for_properties_calculation=True)

pfc.conf_PFC_from_amplitudes()
final_strain = pfc.conf_strain_to_equilibrium()

#Finding elastic constants, starting with mu
pfc_strained = copy.deepcopy(pfc)
f0 = pfc.calc_free_energy()/pfc.volume


strain_limit = 0.05
shear_strains=np.linspace(-strain_limit,strain_limit,100)
free_energies = np.zeros_like(shear_strains)
for n in range(len(shear_strains)):
    shear_strain=shear_strains[n]
    pfc_strained = copy.deepcopy(pfc)
    pfc_strained.conf_apply_distortion(np.array([[0,shear_strain],[0,0.0]]))
    f = pfc_strained.calc_free_energy()/pfc_strained.volume
    free_energies[n] = f-f0

params,_ = curve_fit(lambda x, mu: mu*x**2/2, shear_strains, free_energies)
mu = params[0]
print(f'Eq. mu: {mu:.05f}')
# fig = go.Figure()
# fig.add_trace(go.Scatter(x=shear_strains, y=free_energies, mode='lines', name='Free energy'))
# fig.show
# fig.add_trace(go.Scatter(x=shear_strains, y=mu*shear_strains**2/2, mode='lines', name='Fit'))
plt.plot(shear_strains, free_energies)
plt.plot(shear_strains, mu*shear_strains**2/2, 'r--')
legend = ['Elastic energy density', '$\mu \epsilon^2/2$ (fit)']
plt.legend(legend)
plt.xlabel('Shear strain')
# plt.ylabel('Elastic energy')
plt.grid()

plt.show()


# TO be refined
# mu = 2*(f-f0)/(shear_strain**2)
# print(f'Eq. mu: {mu:.05f}')
# print('Ratio mu eq./mu proto: {:.05f}'.format(mu/pfc.el_mu))


# # gamma
# strain_limit = 0.05
# squeeze_strains = np.linspace(-strain_limit,strain_limit,100)
# free_energies = np.zeros_like(squeeze_strains)
# for n in range(len(squeeze_strains)):
#     squeeze_strain=squeeze_strains[n]
#     pfc_strained = copy.deepcopy(pfc)
#     pfc_strained.conf_apply_distortion(np.array([[squeeze_strain,0],[0,-squeeze_strain]]))
#     f = pfc_strained.calc_free_energy()/pfc_strained.volume
#     free_energies[n] = f-f0

#Basically zero. But how do we account for that? 
# plt.plot(squeeze_strains, free_energies - 2*mu*squeeze_strains**2)
# plt.show()


# a = 0.01
# pfc_strained.conf_apply_distortion(np.array([[a,0],[0,-a]]))
# # free_energy_density = np.mean(pfc_strained.calc_PFC_free_energy_density_and_chemical_potential()[0])
# # f1 = np.mean(free_energy_density)
# f1 = pfc_strained.calc_free_energy()/pfc_strained.volume

# gamma = (f1-f0 - 2*mu*a**2)/(a**2)
# # print('Energy difference: ', f1-f0)
# # print('energy accounted for: ', 2*mu*a**2)
# print(f'Eq. gamma: {gamma:.05f}')
# print('Ratio gamma eq./gamma proto: {:.05f}'.format(gamma/pfc.el_gamma))

# # compression_strain1=0.001
# # pfc_strained.conf_apply_distortion(np.array([[compression_strain1,0],[0,-compression_strain1]]))
# # # print('Volume ratio after compression: {:.05f}'.format(pfc_strained.volume/pfc.volume))
# # f1 = pfc_strained.calc_free_energy()/pfc_strained.volume - f0

# # pfc_strained = copy.deepcopy(pfc)
# # compression_strain2=0.002
# # pfc_strained.conf_apply_distortion(np.array([[compression_strain2,0],[0,-compression_strain2]]))
# # f2 = pfc_strained.calc_free_energy()/pfc_strained.volume - f0

# # # gamma = (f-f0 - 2*mu*compression_strain**2)/(compression_strain**2)
# # gamma = (f2 - f1)/(compression_strain2**2-compression_strain1**2) - 2*mu
# # print(f'Eq. gamma: {gamma:.05f}')
# # print('Ratio gamma eq./gamma proto: {:.05f}'.format(gamma/pfc.el_gamma))

# # Applying a compression strain to find lambda  

# compression_strain = -0.001
# pfc_strained = copy.deepcopy(pfc)
# pfc_strained.conf_apply_distortion(np.array([[compression_strain,0],[0,compression_strain]]))
# f1 = pfc_strained.calc_free_energy()/pfc_strained.volume

# el_lambda = (f1-f0 - 2*mu*compression_strain**2 - gamma*compression_strain**2)/(2*compression_strain**2)

# print(f'Eq. lambda: {el_lambda:.05f}')
# print('Ratio lambda eq./lambda proto: {:.05f}'.format(el_lambda/pfc.el_lambda))