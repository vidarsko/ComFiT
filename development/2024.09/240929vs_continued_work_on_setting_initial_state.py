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

pfc = cf.PhaseFieldCrystal2DTriangular(1,1, for_properties_calculation=False)
# pfc = cf.PhaseFieldCrystal2DSquare(1,1, for_properties_calculation=False)

# pfc.conf_PFC_from_amplitudes()
# final_strain = pfc.conf_strain_to_equilibrium()

# #Finding elastic constants, starting with mu
# pfc_strained = copy.deepcopy(pfc)
# f0 = pfc.calc_free_energy()/pfc.volume

# def elastic_energy(strain, el_lambda, el_mu, el_gamma):
#     exx, exy, eyy = strain
#     return el_lambda/2*(exx+eyy)**2 + el_mu*(exx**2 + 2*exy**2 + eyy**2) + el_gamma/2*(exx**2 + eyy**2)

# strain_magnitudes = np.linspace(-0.01,0.01,21)
# exx = np.array([[a,0,a] for a in strain_magnitudes]).flatten()
# exy = np.array([[0,a,0] for a in strain_magnitudes]).flatten()
# eyy = np.array([[a,0,-a] for a in strain_magnitudes]).flatten()
# # exx = [ 0.01,  0.00,  0.01, -0.01,  0.00, -0.01, 0.02,  0.00,  0.02, -0.02,  0.00, -0.02]
# # exy = [ 0.00,  0.01,  0.00,  0.00, -0.01,  0.00, 0.00,  0.02,  0.00,  0.00, -0.02,  0.00]
# # eyy = [ 0.01,  0.00, -0.01, -0.01,  0.00,  0.01, 0.02,  0.00, -0.02, -0.02,  0.00,  0.02]   

# free_energies = np.zeros_like(exx)
# for n in range(len(exx)):
#     distortion = [[exx[n], exy[n]],[exy[n], eyy[n]]]
#     pfc_strained = copy.deepcopy(pfc)
#     pfc_strained.conf_apply_distortion(distortion)
#     f = pfc_strained.calc_free_energy()/pfc_strained.volume
#     free_energies[n] = f-f0

# params,_ = curve_fit(elastic_energy, (exx, exy, eyy), free_energies)
# el_lambda, el_mu, el_gamma = params
# print(f'Eq. lambda: {el_lambda:.05f}')
# print(f'Eq. mu: {el_mu:.05f}')
# print(f'Eq. gamma: {el_gamma:.05f}')

# def elastic_energy(U, el_lambda, el_mu, el_gamma):
#     return el_lambda/2*np.sum(U**2) + el_mu*np.sum(U**2) + el_gamma*np.sum(U**4) #Something like, but not exactly, this

# strain_limit = 0.01
# shear_strains=np.linspace(-strain_limit,strain_limit,100)
# free_energies = np.zeros_like(shear_strains)
# for n in range(len(shear_strains)):
#     shear_strain=shear_strains[n]
#     pfc_strained = copy.deepcopy(pfc)
#     pfc_strained.conf_apply_distortion(np.array([[0,shear_strain],[0,0.0]]))
#     f = pfc_strained.calc_free_energy()/pfc_strained.volume
#     free_energies[n] = f-f0

# params,_ = curve_fit(lambda x, mu: mu*x**2/2, shear_strains, free_energies)
# mu = params[0]
# print(f'Eq. mu: {mu:.05f}')
# plt.plot(shear_strains, free_energies)
# plt.plot(shear_strains, mu*shear_strains**2/2, 'r--')
# # plt.show()

# # TO be refined
# mu = 2*(f-f0)/(shear_strain**2)
# print(f'Eq. mu: {mu:.05f}')
# print('Ratio mu eq./mu proto: {:.05f}'.format(mu/pfc.el_mu))


# # gamma
# strain_limit = 0.01
# squeeze_strains = np.linspace(-strain_limit,strain_limit,100)
# free_energies = np.zeros_like(squeeze_strains)
# for n in range(len(squeeze_strains)):
#     squeeze_strain=squeeze_strains[n]
#     pfc_strained = copy.deepcopy(pfc)
#     pfc_strained.conf_apply_distortion(np.array([[squeeze_strain,0],[0,-squeeze_strain]]))
#     f = pfc_strained.calc_free_energy()/pfc_strained.volume
#     free_energies[n] = f-f0

# # Basically zero. But how do we account for that? 
# plt.plot(squeeze_strains, free_energies - 2*mu*squeeze_strains**2)
# params,_ = curve_fit(lambda x, gamma: gamma*x**2, squeeze_strains, free_energies - 2*mu*squeeze_strains**2)
# gamma = params[0]
# print(f'Eq. gamma: {gamma:.05f}')
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