import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scipy as sp


pfc = cf.PhaseFieldCrystal3DBodyCenteredCubic(10,10,10)
pfc.conf_PFC_from_amplitudes()
fig = pfc.plot_field(pfc.psi)
fig.show()