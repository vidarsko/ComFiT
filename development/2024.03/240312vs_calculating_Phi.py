import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np

pfc = cf.PhaseFieldCrystal2DTriangular(1,1)
print('Phi Triangular', pfc.Phi)

pfc = cf.PhaseFieldCrystal2DSquare(1,1)
print('Phi Square', pfc.Phi)

pfc = cf.PhaseFieldCrystal3DBodyCenteredCubic(1,1,1)
print('Phi BCC', pfc.Phi)

pfc = cf.PhaseFieldCrystal3DFaceCenteredCubic(1,1,1)
print('Phi FCC', pfc.Phi)

pfc = cf.PhaseFieldCrystal3DSimpleCubic(1,1,1)
print('Phi SC', pfc.Phi)