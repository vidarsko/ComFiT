import comfit as cf
import numpy as np
import matplotlib.pyplot as plt

bec = cf.BEC(2,xRes=101,yRes=101)
bec.set_harmonic_potential(50)
bec.set_initial_condition_Thomas_Fermi()

print(bec.V_ext)



plt.imshow(bec.psi)
plt.colorbar()
plt.show()