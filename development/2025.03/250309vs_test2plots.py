# import comfit as cf
# import matplotlib.pyplot as plt
import numpy as np

# # Set up a 1D quantum system
# # qm = cf.QuantumMechanics(2, xlim=(-50,50), ylim = (-20,20), xRes=1001, yRes=401, dt=0.1)
# # Set up a 1D quantum system
# qm = cf.QuantumMechanics(3, xlim=(-30,30), ylim = (-10,10), zlim=(-10,10), xRes=151, yRes=51, zRes=51, dt=0.1)

# velocity = 1  # Set the velocity of the wavepacket
# potential = 1  # Set the height of the potential barrier

# # Initialize a Gaussian wavepacket at x=5 with velocity=1
# # qm.conf_initial_condition_Gaussian(position=[5,0], width=1, initial_velocity=[1,0])
# qm.conf_initial_condition_Gaussian(position=[5,0,0], width=1, initial_velocity=[velocity,0,0])

# # Add a potential barrier (optional)
# qm.V_ext = 0.7 * (qm.x > 10) * (qm.x < 12)  # Barrier from x=10 to 12

# # height = max((np.max(abs(qm.psi)), np.max(qm.V_ext)))  # Get the maximum value of the wavefunction

# # fig, ax = qm.plot_complex_field(qm.psi)
# # fig, ax = qm.plot_field(qm.V_ext, colorbar=False, opacity = 0.2,
# #                 colormap='bluewhitered', vlim_symmetric=True, 
# #                 xlim=(0,20), ylim = (-10,10), fig=fig, ax=ax)

# fig, ax = qm.plot_complex_field(qm.psi)
# fig, ax = qm.plot_field(qm.V_ext, colorbar=False, 
#                 vlim_symmetric=False, colormap='Reds',
#                 xlim=(0,20), ylim = (-10,10), zlim=(-10,10), fig=fig, ax=ax)

# fig.update_layout(
#     scene_camera=dict(
#         up=dict(x=0, y=0, z=1),       # Sets which direction is up
#         center=dict(x=0, y=0, z=0),    # Sets the center of the camera's focus
#         eye=dict(x=-0.2, y=1.5, z=1.0)  # Sets the camera position
#     )
# )

# qm.show(fig)

for velocity in np.linspace(2, 3, 11):  # Loop over different velocities
    for potential in np.linspace(0.8, 3, 21):  # Loop over different potentials
        print(velocity, potential)