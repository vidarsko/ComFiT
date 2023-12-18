import comfit as cf
import matplotlib.pyplot as plt

pfc = cf.PhaseFieldCrystal3DFaceCenteredCubic(5,5, 5)

pfc.conf_PFC_from_amplitudes()
print(pfc.A, pfc.B)
#pfc.conf_PFC_from_amplitudes()


# for n in range(30):
#     plt.gcf().clear()
#     pfc.evolve_PFC(100)
#     pfc.plot_field(pfc.psi,colormap='viridis')
#     print(pfc.psi.max(),pfc.psi.min())
#     plt.draw()
#     plt.pause(0.01)