import comfit as cf
import matplotlib.pyplot as plt

qm = cf.QM(1)
qm.dt = 0.05
qm.set_initial_condition_gaussian(0.51*qm.xmax,40)
qm.V_ext = qm.set_harmonic_potential(0.0001)


for n in range(200):
    print(n)
    qm.evolve_schrodinger(200)
    plt.cla()
    qm.plot()
    plt.draw()
    plt.pause(0.01)