import comfit as cf

axis_limit = 3
qm = cf.QuantumMechanics(3,xlim=[-axis_limit,axis_limit],ylim=[-axis_limit,axis_limit],zlim=[-axis_limit,axis_limit],
                            xRes=50,yRes=50,zRes=50)
qm.conf_hydrogen_state(1, 0, 0)
# fig, ax = qm.plot_complex_field(qm.psi)
# fig

axis_limit = 10
qm = cf.QuantumMechanics(3,xlim=[-axis_limit,axis_limit],ylim=[-axis_limit,axis_limit],zlim=[-axis_limit,axis_limit])
qm.conf_hydrogen_state(2, 1, 1)
# fig, ax = qm.plot_complex_field(qm.psi,phase_blob_threshold=0.2)
# fig

for n in range(50):
    qm.evolve_schrodinger(5)
    fig, ax = qm.plot_complex_field(qm.psi,phase_blob_threshold=0.2)
    qm.plot_save(fig, n)

cf.tool_make_animation_gif(n,name="hydrogen_evolution_free")