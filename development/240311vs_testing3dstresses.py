import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import matplotlib.pyplot as plt

pfc = cf.PhaseFieldCrystal3DBodyCenteredCubic(13, 13, 13)
eta = pfc.calc_amplitudes_with_dislocation_ring()
pfc.conf_PFC_from_amplitudes(eta)
pfc.evolve_PFC(10)

stress = pfc.calc_stress_tensor()

pfc.plot_field_in_plane(stress[0]/pfc.el_mu)
plt.show()

# fig = plt.figure()
# axs = fig.subplots(3,3,subplot_kw={'projection': '3d'})

# for i in range(3):
#     for j in range(3):
#         pfc.plot_field_in_plane(pfc.get_sym(stress,i,j)/pfc.el_mu,ax=axs[i,j])