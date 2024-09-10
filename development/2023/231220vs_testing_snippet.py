import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import matplotlib.pyplot as plt

pfc = cf.BoseEinsteinCondensate(3)
pfc.conf_insert_vortex_ring() 
pfc.evolve_relax(100)
pfc.plot_field(abs(pfc.psi))

plt.show()