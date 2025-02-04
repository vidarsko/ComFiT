
import comfit as cf


pfc = cf.PhaseFieldCrystal2DSquare(40,40)
eta = pfc.calc_amplitudes_with_dislocation_dipole()
pfc.conf_PFC_from_amplitudes(eta)
pfc.evolve_PFC(100)

pfc.plot_lib = 'plotly'

fig = pfc.plot_field(pfc.psi)
pfc.plot_save(1, fig)