import comfit as cf
import matplotlib.pyplot as plt

pfc = cf.PFCtri(13,7)

eta = pfc.conf_insert_dislocation_dipole(dislocation_type=1)
pfc.calc_from_amplitudes(eta)

#ax = plt.gcf().add_subplot(111)

for n in range(10):
    plt.gcf().clear()
    pfc.plot_field(pfc.psi)
    pfc.evolve_PFC(200)
    plt.draw()
    plt.pause(0.05)
