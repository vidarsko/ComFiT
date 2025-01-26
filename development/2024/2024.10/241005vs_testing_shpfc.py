import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import comfit as cf
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scipy as sp

# # import comfit as cf
# # from activePFC import activePFC
# import matplotlib.pyplot as plt
# import numpy as np


for activity in np.linspace(0,0,1):
    pfc = cf.PhaseFieldCrystal2DSquare(20, 20, plot_lib='matplotlib')
    # pfc.conf_PFC_from_amplitudes()
    eta = pfc.calc_amplitudes_with_dislocation_dipole()
    pfc.conf_PFC_from_amplitudes(eta)
    pfc.evolve_PFC(500)
    # print(pfc.plot_lib)

    for n in range(200):
        plt.gcf().clear()
        pfc.evolve_PFC_hydrodynamic(20,gamma_S = 2**-6, rho0 = 2**-6)
        pfc.plot_PFC()
        cf.tool_save_plot_matplotlib(n)
        # plt.show()
        # plt.pause(0.05)

    cf.tool_make_animation_gif(n)