import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit as cf
import matplotlib.pyplot as plt

qm = cf.QuantumMechanics(3, xRes=51, yRes=51,zRes=51)
qm.conf_initial_condition_gaussian(position=qm.rmid,width=5)
qm.V_ext = qm.conf_harmonic_potential(0.0001)
qm.evolve_schrodinger(100)
qm.plot_complex_field(qm.psi)
plt.show()
# print(qm.psi)

# for i in range(150):
#     # qm.plot_field(abs(qm.psi))
#     qm.plot_complex_field(qm.psi)
#     qm.evolve_schrodinger(100)
#     plt.draw()
#     plt.pause(0.01)
    # cf.tool_save_plot(i)
# cf.tool_make_animation(i)

