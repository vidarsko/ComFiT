import comfit as cf
import matplotlib.pyplot as plt


bec = cf.BEC(3,xRes=21,yRes=21,zRes=21)
bec.set_initial_condition_disordered()
bec.evolve_relax_BEC(200)

bec.plot_field(abs(bec.psi))
plt.show()