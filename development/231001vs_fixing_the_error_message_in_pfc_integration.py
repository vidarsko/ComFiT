import comfit as cf
import matplotlib.pyplot as plt

pfc = cf.PhaseFieldCrystal2DTriangular(7,3)
omega_f = pfc.r + (1 - pfc.calc_k2())**2
integrating_factors_f = pfc.calc_evolution_integrating_factors_f()


print(integrating_factors_f[2][0,0])
pfc.plot_fourier_field(integrating_factors_f[2])
#pfc.plot_fourier_field(omega_f)
plt.show()