import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import plotly.graph_objects as go

# pfc = cf.PhaseFieldCrystal2DTriangular(12,5)
# pfc = cf.PhaseFieldCrystal1DPeriodic(3,micro_resolution=[15])
# pfc = cf.PhaseFieldCrystal3DBodyCenteredCubic(3,3,3)
pfc.plot_lib = 'plotly'

pfc.conf_PFC_from_amplitudes()
fig = pfc.plot_field(pfc.psi)
fig.show()