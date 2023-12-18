import comfit as cf
import matplotlib.pyplot as plt

pfc = cf.PhaseFieldCrystal2DTriangular(20,10)
print(pfc.a0)
print(pfc.q)

eta = pfc.calc_amplitudes_with_dislocation()
pfc.conf_PFC_from_amplitudes(eta)

#fig = plt.figure()
#ax = plt.gcf().add_subplot(111)

#We are in business
pfc.plot_field(pfc.psi)
#plt.colorbar(plt.gci())
plt.show()