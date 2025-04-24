import comfit as cf

qm = cf.QuantumMechanics(1, xlim=[-100,100], xRes=201, dt=0.1)


alpha=10
V0=1.0
qm.V_ext = (qm.x > 0)*(qm.x<alpha) * V0


V = qm.calc_coarse_grain(qm.V_ext,1)

fig, ax = qm.plot_field(V, title='Coarse-grained potential')
qm.show(fig)
