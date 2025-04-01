import comfit as cf

pfc = cf.PhaseFieldCrystal2DTriangular(4,4)
pfc.conf_PFC_from_amplitudes()
fig, ax = pfc.plot_field(pfc.psi)
pfc.show(fig)