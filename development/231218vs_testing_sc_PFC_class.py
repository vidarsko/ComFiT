import comfit as cf
import matplotlib.pyplot as plt

#pfc = cf.PhaseFieldCrystal3DFaceCenteredCubic(5,5, 5)
pfc = cf.PhaseFieldCrystal3DSimpleCubic(31,31, 31)
eta = pfc.calc_amplitudes_with_dislocation_ring(normal_vector=[0,0,1])
pfc.conf_PFC_from_amplitudes(eta)
#print(pfc.A, pfc.B, pfc.C)
#print(pfc.el_lambda, pfc.el_mu, pfc.el_gamma, pfc.el_nu)
#pfc.conf_PFC_from_amplitudes()

pfc.evolve_PFC(50)
#pfc.plot_field(pfc.psi,colormap='viridis',number_of_layers=1)
#print(pfc.psi)
pfc.plot_field_in_plane(pfc.psi,normal_vector=[0,1,0])
plt.show()

# for n in range(30):
#     plt.gcf().clear()
#     pfc.evolve_PFC(100)
#     pfc.plot_field(pfc.psi,colormap='viridis')
#     print(pfc.psi.max(),pfc.psi.min())
#     plt.draw()
#     plt.pause(0.01)