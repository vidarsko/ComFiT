import comfit as cf     # type: ignore


bec = cf.BoseEinsteinCondensate(3)
bec.conf_insert_vortex_ring()
bec.evolve_relax(200)

fig1 = bec.plot_complex_field(bec.psi, phase_blob_threshold=0.9)
fig2 = bec.plot_field(abs(bec.psi), number_of_layers=3)
fig = cf.tool_make_subplots(1, 2, fig1, fig2)
fig.show()
