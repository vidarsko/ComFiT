

import comfit as cf

qm = cf.QuantumMechanics(1)
qm.conf_initial_condition_Gaussian()

# interval = qm.calc_region_interval(5,10)
# print(qm.calc_integrate_field(abs(qm.psi)**2))

def integrate_interval():
    


interval = (qm.x>5) & (qm.x<10)
np.sum(interval*qm.dx*abs(qm.psi)**2)