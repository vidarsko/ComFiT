import comfit as cf
import matplotlib.pyplot as plt

pfc = cf.PFCtri(20,20)
print(pfc.a0)
print(pfc.q)

pfc.calc_from_amplitudes()

#We are in business
plt.imshow(pfc.Psi)
plt.show()