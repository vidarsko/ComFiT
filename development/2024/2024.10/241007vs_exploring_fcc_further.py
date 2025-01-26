import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scipy as sp


pfc = cf.PhaseFieldCrystal3DFaceCenteredCubic(1,1,1)
pfc.conf_PFC_from_amplitudes()
pfc.evolve_PFC(100)

