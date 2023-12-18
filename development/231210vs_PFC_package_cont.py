import comfit as cf
import matplotlib.pyplot as plt

pfc = cf.PhaseFieldCrystal2DTriangular(21,14, dt=0.1)
eta = pfc.calc_amplitudes_with_dislocation_dipole()
pfc.conf_PFC_from_amplitudes(eta)
pfc.evolve_PFC(200)
pfc.plot_field(pfc.psi)
plt.show()
