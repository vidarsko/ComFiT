import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import matplotlib.pyplot as plt


bec = cf.BoseEinsteinCondensate(3,xRes=31,yRes=31,zRes=31)
bec.conf_initial_condition_disordered()
bec.evolve_relax(100)

bec.plot_complex_field(bec.psi, plot_method='phase_blob')

plt.show()