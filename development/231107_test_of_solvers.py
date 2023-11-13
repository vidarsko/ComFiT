import comfit as cf
import numpy as np
import matplotlib.pyplot as plt

bec = cf.BEC(2,xRes=3,yRes=3)

bec.psi = np.arange(9).reshape((3,3))
bec.psi_f =np.fft.fft2(bec.psi)

integrating_gactors = bec.calc_evolution_integrating_factors_dGPE_f()

new_t,new_tf =bec.evolve_ETD2RK_loop_test(integrating_gactors,bec.calc_nonlinear_evolution_term_f,bec.psi,bec.psi_f)
new,new_f = bec.evolve_ETD2RK_loop(integrating_gactors,bec.calc_nonlinear_evolution_term_f,bec.psi,bec.psi_f)


#print(new)
#print(new_t)