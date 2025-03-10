import comfit as cf
import matplotlib.pyplot as plt
import numpy as np

# # Set up a 1D quantum system
qm = cf.QuantumMechanics(3, xlim=(-30,30), ylim = (-10,10), zlim=(-10,10), xRes=151, yRes=51, zRes=51, dt=0.05)

velocity = 2.4
for potential in np.linspace(0.8, 3, 21):  # Loop over different potentials
    phase_blob_threshold = 0.5

    # Initialize a Gaussian wavepacket at x=5 with velocity=1
    qm.conf_initial_condition_Gaussian(position=[5,0,0], width=1, initial_velocity=[velocity,0,0])

    # Add a potential barrier (optional)
    qm.V_ext = potential * (qm.x > 10) * (qm.x < 12)  # Barrier from x=10 to 12

    # height = max((np.max(abs(qm.psi)), np.max(qm.V_ext)))  # Get the maximum value of the wavefunction
    # Optional: Animate it
    for n in range(305):
        print(n)
        fig, ax = qm.plot_complex_field(qm.psi, phase_blob_threshold=phase_blob_threshold)
        fig, ax = qm.plot_field(qm.V_ext, colorbar=False, 
                    vlim_symmetric=False, colormap='Reds',
                    xlim=(0,20), ylim = (-10,10), zlim=(-10,10), fig=fig, ax=ax)

        fig.update_layout(
        scene_camera=dict(
            up=dict(x=0, y=0, z=1),       # Sets which direction is up
            center=dict(x=0, y=0, z=0),    # Sets the center of the camera's focus
            eye=dict(x=-0.2, y=1.5, z=1.0)  # Sets the camera position
        )
        )
        qm.plot_save(fig, n, ID=f'{velocity}_{potential}_{phase_blob_threshold}')  # Save the figure

        qm.evolve_schrodinger(1)
    cf.tool_make_animation_movie(n, ID=f'{velocity}_{potential}_{phase_blob_threshold}')  # Creates animation

# cf.tool_make_animation_movie(44, ID=f'{velocity}_{potential}')  # Creates animation