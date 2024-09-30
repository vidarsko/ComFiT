import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scipy as sp
from comfit.plot.plot_field_matplotlib import plot_field_matplotlib


pfc = cf.PhaseFieldCrystal2DSquare(3,3)
pfc.conf_PFC_from_amplitudes()
pfc.conf_strain_to_equilibrium()
pfc.conf_apply_distortion(np.array([[0,-0.05],[0,0.0]]))
plot_field_matplotlib(pfc,pfc.psi, X=pfc.X, Y=pfc.Y)
plt.axis('tight')
plt.show()