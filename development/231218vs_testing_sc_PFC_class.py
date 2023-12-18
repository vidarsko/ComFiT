import comfit as cf
import matplotlib.pyplot as plt

pfc = cf.PhaseFieldCrystal3DSimpleCubic(1,1, 1)
#pfc = cf.PhaseFieldCrystal3DSimpleCubic(1,1, 1,psi0=-0.225,r=-0.2)

pfc.conf_PFC_from_amplitudes()
print(pfc.A, pfc.B, pfc.C)
#pfc.conf_PFC_from_amplitudes()


# for n in range(30):
#     plt.gcf().clear()
#     pfc.evolve_PFC(100)
#     pfc.plot_field(pfc.psi,colormap='viridis')
#     print(pfc.psi.max(),pfc.psi.min())
#     plt.draw()
#     plt.pause(0.01)