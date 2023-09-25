import comfit as cf
import matplotlib.pyplot as plt

pfc = cf.PFCtri(20,10)
print(pfc.a0)
print(pfc.q)

eta = pfc.conf_insert_dislocation()
pfc.calc_from_amplitudes(eta)

ax = plt.subplot()

#We are in business
pfc.plot_field(pfc.psi,ax)
#plt.colorbar(plt.gci())
plt.show()