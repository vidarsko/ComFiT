import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scipy as sp

# pfc = cf.PhaseFieldCrystal2DTriangular(1,1, type_of_evolution='unconserved')

# from unittest import TestCase


pfc = cf.PhaseFieldCrystal2DTriangular(31, 16)
x1=pfc.xmax/3 
y1=pfc.ymax/2
x2=2*pfc.xmax/3
y2=pfc.ymax/2

eta = pfc.calc_amplitudes_with_dislocation_dipole(
    dislocation_type=1,
    x1=x1, y1=y1,
    x2=x2, y2=y2)
pfc.conf_PFC_from_amplitudes(eta)
pfc.evolve_PFC(100)

# Check if there are two dislocations
dislocation_nodes = pfc.calc_dislocation_nodes()
assert len(dislocation_nodes) == 2

pfc.evolve_PFC_mechanical_equilibrium(350,Delta_t=1)

# Check that the dislocation\s have annihilated
dislocation_nodes = pfc.calc_dislocation_nodes()
assert len(dislocation_nodes) == 0
# TestCase.assertEqual(len(dislocation_nodes),0)