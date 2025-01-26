import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

bec = cf.BoseEinsteinCondensate(3,xRes=65,yRes=65,zRes=65,gamma=0.005,dt=0.1)

bec.conf_insert_vortex_ring(radius=10, normal_vector=[0, 1, 0])

bec.conf_insert_vortex_ring(radius=20, normal_vector=[0, -1, 0])

bec.evolve_relax(20)

N = 100
for n in range(N):
    bec.evolve_dGPE(30)
    fig,ax=bec.plot_field(abs(bec.psi)**2, vlim=[0,1])
    
    cf.tool_save_plot_matplotlib(n)


cf.tool_make_animation_gif(n)