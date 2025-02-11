import comfit as cf     # type: ignore


bec = cf.BoseEinsteinCondensate(3,xRes=51,yRes=51,zRes=51)
bec.conf_initial_condition_disordered()

bec.evolve_relax(100)
# fig1 = bec.plot_complex_field(bec.psi)
# fig2 = bec.plot_field(abs(bec.psi))
# fig = cf.tool_make_subplots(1, 2, fig1, fig2)
# bec.evolve_dGPE(2000)

# fig3 = bec.plot_complex_field(bec.psi)
# fig4 = bec.plot_field(abs(bec.psi))
# fig = cf.tool_make_subplots(2, 2, fig1, fig2, fig3, fig4)
# fig.show()

for i in range(100):
    bec.evolve_dGPE(100)
    fig1 = bec.plot_complex_field(bec.psi)
    fig2 = bec.plot_field(abs(bec.psi))
    fig = cf.tool_make_subplots(1, 2, fig1, fig2)
    bec.plot_save(i, fig)
cf.tool_make_animation_gif(i)
