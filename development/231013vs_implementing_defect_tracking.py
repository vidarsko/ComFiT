import comfit as cf

import comfit as cf
import matplotlib.pyplot as plt


bec = cf.BEC(2,xRes=61,yRes=61)
bec.set_initial_condition_disordered()
bec.evolve_relax_BEC(200)

