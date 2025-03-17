import comfit as cf
import matplotlib.pyplot as plt
import numpy as np

# Set up a 1D quantum system
qm = cf.QuantumMechanics(1, xlim=[-50,50], xRes=1001, dt=0.1)

# Initialize a Gaussian wavepacket at x=5 with velocity=1
qm.conf_initial_condition_Gaussian(position=5, width=1, initial_velocity=1)

# Add a potential barrier (optional)
qm.V_ext = 0.5 * (qm.x > 10) * (qm.x < 12)  # Barrier from x=10 to 12

height = np.max(abs(qm.psi))  # Get the maximum value of the wavefunction

# Optional: Animate it
for n in range(61):
    qm.evolve_schrodinger(5)
    fig, ax = qm.plot_complex_field(qm.psi)
    qm.plot_field(qm.V_ext, fig=fig, ax=ax, ylim=[0,height], xlim=[0,20])
    qm.plot_save(fig, n)
cf.tool_make_animation_gif(n)  # Creates animation.gif