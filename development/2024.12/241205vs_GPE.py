import comfit as cf

bec = cf.BoseEinsteinCondensate(2)
bec.conf_initial_condition_Thomas_Fermi()

for n in range(100):
    bec.evolve_dGPE(1)
    fig1 = bec.plot_complex_field(bec.psi, colorbar=True)
    fig2 = bec.plot_field(abs(bec.psi), vlim=[0,0.8221])

    fig = cf.tool_make_subplots(1,2,fig1,fig2)
    fig.update_layout(
    autosize=False,
    width=1000,
    height=400,
    )
    cf.tool_save_plot(n, fig)

cf.tool_make_animation_gif(n)