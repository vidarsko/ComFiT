import comfit as cf
import matplotlib.pyplot as plt

pfc = cf.PFCtri(20,10)

eta = pfc.conf_insert_dislocation(dislocation_type=1)
pfc.calc_from_amplitudes(eta)

pfc.plot_field(pfc.psi)
plt.show()