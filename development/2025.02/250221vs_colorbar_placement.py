import comfit as cf
import numpy as np
import plotly.graph_objects as go


# dt = 0.1
# qm = cf.QuantumMechanics(2, xlim=[-15,15], ylim=[-15,15], xRes=201, yRes=201, dt=dt)

# initial_pos = [5,5]
# k = 0.01
# qm.V_ext = k*(qm.x**2 + qm.y**2)
# qm.conf_initial_condition_Gaussian(position=initial_pos, width=1)

# def plot_qm_vs_classical(T):
#     dt = 0.1 
#     qm = cf.QuantumMechanics(2, xlim=[-30,30], ylim=[-30,30], xRes=201, yRes=201, dt=dt)

#     initial_pos = [5,5]
#     k = 0.01
#     V_ext = 1/2 * k*(qm.x**2 + qm.y**2)

#     qm.V_ext = V_ext
#     qm.conf_initial_condition_Gaussian(position=initial_pos, width=1)
#     # qm.dt = 0.1

#     timesteps = int(T / dt)

#     qm.evolve_schrodinger(timesteps)
#     fig, ax = qm.plot_complex_field(qm.psi)

#     pos = np.array(initial_pos)
#     vel = np.array([0, 0])
#     for _ in range(timesteps):
#         a = - k * pos
#         vel = vel + a * dt
#         pos = pos + vel * dt
        
#     fig.add_trace(go.Scatter(x=[pos[0]], y=[pos[1]], mode='markers', name='Classical solution'))

#     fig.show()


# plot_qm_vs_classical(0)
# Initialize a 1-Dimensional system 
dt = 0.1
qm = cf.QuantumMechanics(1,xlim=[-10,10], dt=dt) 

# Sets the initial condition to a Gaussian wavepacket centered at x=0 with width 1
qm.conf_initial_condition_Gaussian(position=0, width=1) 
psi = qm.psi # Get the wavefunction at time t=0

fig, ax = qm.plot_complex_field(psi)
fig.show()