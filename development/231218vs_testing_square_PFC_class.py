import comfit as cf
import matplotlib.pyplot as plt

pfc = cf.PhaseFieldCrystal2DSquare(21,21)
eta = pfc.calc_amplitudes_with_dislocation_dipole()
pfc.conf_PFC_from_amplitudes(eta)
#pfc.conf_PFC_from_amplitudes()

for n in range(30):
    plt.gcf().clear()
    pfc.evolve_PFC(100)
    pfc.plot_field(pfc.psi,colormap='viridis')
    print(pfc.psi.max(),pfc.psi.min())
    plt.draw()
    plt.pause(0.01)