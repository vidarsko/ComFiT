import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scipy as sp

pfc = cf.PhaseFieldCrystal2DSquare(10,10)
pfc.conf_PFC_from_amplitudes()

fig1 = pfc.plot_field(pfc.psi, axis_equal=True)

pfc2 = cf.PhaseFieldCrystal3DBodyCenteredCubic(2,2,2)
pfc2.conf_PFC_from_amplitudes()

fig2 = pfc2.plot_field(pfc2.psi)

fig = cf.tool_make_subplots(2,2,fig2,fig1,fig1,fig1)
fig.show()