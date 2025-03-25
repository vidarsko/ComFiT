import numpy as np
import matplotlib.pyplot as plt
import comfit as cf

# Initialize a 1-Dimensional system 
dt = 0.1
xRes = 201
xlim = [-10,10]
dx = (xlim[1]-xlim[0])/(xRes-1)

qm = cf.QuantumMechanics(1, xlim=xlim, dt=dt, xRes=xRes)
qm.conf_initial_condition_Gaussian(position=0, width=1, initial_velocity=1)
psi_initial = qm.psi


def integrate_interval(field: np.ndarray, x_start: float, x_end: float, qm: cf.QuantumMechanics) -> float:
    """Integrates a field over a given interval

    Args:
        field (np.ndarray): a field to integrate
        x_start (float): start of the interval
        x_end (float): end of the interval
        qm (cf.QuantumMechanics): quantum mechanics object

    Returns:
        float: the integral of the field over the interval
    """
    x = qm.x
    dx = x[1] - x[0]
    is_in_interval = (x >= x_start) & (x <= x_end) # Boolean array with True where x is in the interval
    return np.sum(field[is_in_interval]) * dx


# Initialize a 1-Dimensional system 
dt = 0.1
xRes = 1001
xlim = [-50,50]
dx = (xlim[1]-xlim[0])/(xRes-1)

qm = cf.QuantumMechanics(1, xlim=xlim, dt=dt, xRes=xRes)
qm.conf_initial_condition_Gaussian(position=-10, width=1, initial_velocity= 1)
# qm.V_ext = 1
# V = 1 * (qm.x > 0 ) * (qm.x < 5)
# qm.V_ext = V
fig, ax = qm.plot_field(qm.V_ext, size=(400, 400))
qm.show(fig)