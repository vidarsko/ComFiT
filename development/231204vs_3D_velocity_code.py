import comfit as cf
import matplotlib.pyplot as plt

bec = cf.BoseEinsteinCondensate(3,xRes=31,yRes=31,zRes=31)
bec.conf_insert_vortex_ring()
bec.evolve_relax(300)

psi_old = bec.psi

bec.evolve_relax(1)
dt_psi = (bec.psi - psi_old)/bec.dt

psi = [np.real(bec.psi),np.imag(bec.psi)]
dt_psi = [np.real(dt_psi),np.imag(dt_psi)]

velocity_field = bec.calc_defect_velocity_field(dt_psi)
