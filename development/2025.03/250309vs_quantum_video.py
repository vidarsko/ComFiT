import comfit as cf
import matplotlib.pyplot as plt
import numpy as np

# Set up a 1D quantum system
qm = cf.QuantumMechanics(2, xlim=(-50,50), ylim = (-20,20), xRes=1001, yRes=401, dt=0.1)

velocity = 1.5  # Set the velocity of the wavepacket
potential = 3  # Set the height of the potential barrier

# Initialize a Gaussian wavepacket at x=5 with velocity=1
qm.conf_initial_condition_Gaussian(position=[5,0], width=1, initial_velocity=[velocity,0])

# Add a potential barrier (optional)
qm.V_ext = potential * (qm.x > 10) * (qm.x < 12)  # Barrier from x=10 to 12

# height = max((np.max(abs(qm.psi)), np.max(qm.V_ext)))  # Get the maximum value of the wavefunction
# Optional: Animate it
for n in range(305):
    print(n)
    qm.evolve_schrodinger(1)
    fig, ax = qm.plot_complex_field(qm.psi)
    fig, ax = qm.plot_field(qm.V_ext, colorbar=False, opacity = 0.2,
                    colormap='bluewhitered', vlim_symmetric=True, 
                    xlim=(0,20), ylim = (-10,10), fig=fig, ax=ax)
    
    qm.plot_save(fig, n, ID=f'{velocity}_{potential}')  # Save the figure
cf.tool_make_animation_movie(n, ID=f'{velocity}_{potential}')  # Creates animation