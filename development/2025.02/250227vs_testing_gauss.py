import comfit as cf
import numpy as np

# Create quantum mechanics instance
qm = cf.QuantumMechanics(1)  # 1D quantum system

# Set up initial Gaussian wavepacket
qm.conf_initial_condition_Gaussian(
    position=0.0,      # Center position
    width=2,         # Width of Gaussian
    initial_velocity=2.0  # Initial momentum
)

# Create animation
for n in range(100):
    # Evolve quantum system
    qm.evolve_schrodinger(1)  # Single time step evolution

    # Plot and save frame - showing probability density
    fig, ax = qm.plot_field(np.abs(qm.psi)**2, title=f'Time: {qm.time:.2f}')
    qm.plot_save(fig, n)

# Generate the final animation from saved frames
cf.tool_make_animation_gif(n)  # Number of frames