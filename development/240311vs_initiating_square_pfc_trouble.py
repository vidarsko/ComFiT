import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf

pfc = cf.PhaseFieldCrystal2DSquare(30,30)
print(pfc.x.shape)
print(pfc.xRes)