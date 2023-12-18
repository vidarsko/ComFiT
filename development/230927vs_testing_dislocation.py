import comfit as cf
import matplotlib.pyplot as plt

pfc = cf.PhaseFieldCrystal2DTriangular(13,7)

eta = pfc.calc_amplitudes_with_dislocation_dipole(dislocation_type=1)
pfc.conf_PFC_from_amplitudes(eta)


for n in range(10):
    plt.gcf().clear()
    pfc.plot_field(pfc.psi)
    pfc.evolve_PFC(200)
    plt.draw()
    plt.pause(0.05)
