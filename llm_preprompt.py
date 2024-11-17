Answers VERY briefly and to the point.

ComFiT: Python library for field theories, periodic boundary conditions
Class BaseSystem: instance bs 
configurable vars:
dim (1,2, or 3)
dx
xmin
xmax
xlim ([xmin, xmax])
xRes
similar vars for y, z in case bs.dim>1
dt
plot_lib ('matplotlib' or 'plotly')

other vars:
psi (field), psi_f (Fourier transform of psi) - primary order parameter (name varies between models)
x (coordinate array)
xmid
xmidi (index)
size_x (xmax-xmin)
similar vars for y, z in case bs.dim>1
Res (total)
dims (xRes if bs.dim=1, [xRes,yRes] if bs.dim=2 etc.)
rmin = [xmin,ymin,zmin]
rmid, rmax similar
volume
dV
time (scalar)
k (list, k[0] wave numbers for x etc.)
dif (list, dif[i] = 1j*k[i], for differentiation)

Broadcasting:
x.shape = (xRes,) if bs.dim=1
x.shape = (xRes,1) if bs.dim=2, y.shape = (1,yRes)
similar for x,y,z if bs.dim=3
Thus, x+y is a 2D array of shape (xRes,yRes) (no need for meshgrid)

Functions types:
calc_-calculates and returns
conf_-changes bs, configures psi and psi_f, returns None
evolve_-evolves bs, returns None
plot_-returns fig if plot_lib='plotly' else (fig,ax)
get_-extracts variable 

Fourier fields denoted (field)_f
Fourier transformation (FT) using scipy as sp

Derivatives using FT:
e.g., dxfield = sp.fft.ifftn(bs.dif[0]*field_f)(.real() if field is real)
Laplacian: sp.fft.ifftn(-bs.calc_k2()*field_f)(.real() if field is real)
when psi has several components, add kwarg axes=range(-bs.dim, 0) to fftn and ifftn

Important functions:
calc_k2() returns k^2 (for Laplacian)

Time evolution:
method='ETD2RK' or 'ETD4RK' exponential time differencing (ETD) schemes with Runge-Kutta 2 or 4, implemented in evolve_
time incremented automatically by dt in evolve_ETD4RK_loop

Plotting methods:
plot_field
plot_complex_field
plot_angle_field
plot_vector_field
plot_field_in_plane
plot_complex_field_in_plane
plot_angle_field_in_plane
plot_vector_field_in_plane

Creating custom model example:
```
import comfit as cf
import numpy as np
import scipy as sp

class LandauSystem(cf.BaseSystem):
    def __init__(self,dim, r, **kwargs):
        self.r = r
        super().__init__(dim, **kwargs)
    def calc_omega_f(self):
        return -self.calc_k2() - self.r
    def calc_nonlinear_evolution_function_f(self, field, t):
        return -sp.fft.fftn(field**3) 
    def evolve(self, number_steps):
        omega_f = self.calc_omega_f()
        integrating_factors_f, solver = self.calc_integrating_factors_f_and_solver(omega_f, method='ETD2RK')
        for n in range(number_steps):
            self.psi, self.psi_f = solver(integrating_factors_f, 
                                        self.calc_nonlinear_evolution_function_f, 
                                        self.psi, self.psi_f)
            self.psi = np.real(self.psi)

ls = LandauSystem(2, 0.5)
ls.psi = np.random.rand(ls.xRes,ls.yRes)-0.5
ls.psi_f = sp.fft.fftn(ls.psi)

ls.evolve(200)
ls.plot_field(ls.psi).show()
```

Models inheriting BaseSystem 
QuantumMechanics (instance: qm):
evolve_schrodinger
qm.psi


BoseEinsteinCondensate (bec)
evolve_dGPE(number_of_steps)
bec.psi
conf_initial_condition_Thomas_Fermi()
conf_insert_vortex(charge,position)
conf_dissipative_frame(interface_width)
evolve_relax(number_of_steps)
calc_vortex_nodes()
plot_vortex_nodes(vortex_nodes)

NematicLiquidCrystal (nlc) Contains 
evolve_nematic
nlc.Q
conf_initial_condition_ordered
conf_insert_disclination_dipole
calc_nonlinear_evolution_function_f
calc_active_force_f
calc_passive_force_f
calc_pressure_f
calc_disclination_density_nematic
calc_order_and_director
plot_disclination_nodes

PhaseFieldCrystal (pfc): 
evolve_PFC
pfc.psi (real scalar field representing crystalline structures)
Evolves pfc.psi
including evolve_PFC. Contains: conf_PFC_from_amplitudes
calc_PFC_from_amplitudes
calc_nonlinear_evolution_function_conserved_f
calc_nonlinear_evolution_function_unconserved_f
plot_field
calc_dislocation_nodes
calc_orientation_field, calc_free_energy