import comfit as cf
import numpy as np
import matplotlib.pyplot as plt

bec = cf.BEC(2,xRes=101,yRes=101)
bec.set_harmonic_potential(50)
bec.set_initial_condition_Thomas_Fermi()

integrating_factors = bec.calc_evolution_integrating_factors_dGPE_f()
beta =bec.calc_nonlinear_evolution_term_dGPE_f(bec.psi)



plt.imshow(np.abs(integrating_factors[2]))
plt.colorbar()
plt.show()