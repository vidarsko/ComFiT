import comfit as cf
import matplotlib.pyplot as plt

pfc = cf.PFCtri(20,10)
print(pfc.a0)
print(pfc.q)

pfc.calc_from_amplitudes()

ax = plt.subplot()

#We are in business
pfc.plot_field(pfc.psi,ax)
plt.show()