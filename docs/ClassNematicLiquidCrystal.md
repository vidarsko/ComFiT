# Active Nematics

## Variables

The primary variables are the symmetric traceless tensor $Q$ and the
velocity field $\\vec{u}$

``` {.python language="Python"}
nem.Q 
    nem.u
```

The $Q$ tensor is given by the nematic director as 

$$
Q_{ij} = S (n_i n_j - \frac{1}{d} \delta_{ij})
$$ 

in $d$ dimensions. For
the current version only $d=2$ is implemented. To take advantage of its
symmetric nature we have saved $Q$ as a vector field, which in two
dimensions takes the form $$\\vec{Q} =[ Q_{xx},Q_{xy}. ]$$ We can
translate from tensor indexes to the right value stored in the vector by
using the function

``` {.python language="Python"}
get_sym(self,Q,i,j)
```

This returns the element $Q_{ij}$ of a symmetric and traceless tensor
field. In addition to this we also have the function

``` {.python language="Python"}
get_anti_sym(self,omega,i,j)
```

so that we can optimally store the antisymetric tensors as well. In two
dimensions these only have one independent component, which is stored as
a scalar field.

## Model

We model the active nematic using a set of coupled differential
equations, namely the Edvard-Beris equation coupled to the Stokes
equation
[@marchetti2013hydrodynamics; @genkin2017topological; @nejad2020memory; @angheluta2021role]

$$
\begin{aligned}
\partial_t Q + \mathbf u\cdot \nabla Q +\Omega Q -Q \Omega &=-\gamma^{-1}H,
\\
(\Gamma- \eta \nabla^2 )\mathbf u &= -\nabla P + \nabla \cdot \sigma^a(Q) + \nabla \cdot \sigma^p, \\
    \nabla \cdot \mathbf u &= 0.
\end{aligned}
$$ 

Here $2\Omega_{ij} = \partial_i u_j - \partial_j u_i$ is
the vorticity tensor, $P$ is the pressure and we have the active stress
$\sigma^a = \alpha Q$. The vorticity tensor is calculated by the
function

``` {.python language="Python"}
calc_vorticity_tensor(self)
```

Note that the velocity has to be updated before this function is called.
The calculation of the pressure and velocity is described in section
[7.4](#sec:nem_vel){reference-type="ref" reference="sec:nem_vel"}. Since
the active stress is proportional to $Q$ we have not included any
function to calculate it.

The molecular field $H$ is given as
$$H_{ij} =  \frac{\delta \mathcal{F}}{\delta Q_{ij}} =- K \nabla^2 Q_{ij} -A(B - 2Q^2_{kk})Q_{ij}$$
in two dimensions. We have here used the free energy
$$\mathcal F = \int \left( K |\nabla Q|^2 - \frac{A}{2} \left[ B \text{Tr}(Q^2) -\text{Tr}(Q^2)^2   \right] \right),$$
where it is assumed that there is a single Frank elastic constant $K$.
For the passive stress we use
$$\sigma_{ij} = Q_{ik}H_{kj}-H_{ik}Q_{kj} - K (\partial_i Q_{kl})(\partial_jQ_{kl}).$$
Note that we use the convention that the dot product between a vector
and a tensor contracts the last component when calculating the
divergence of this. The first two terms is the asymmetric stress, while
the second term is the Ericksen stress. Terms due to flow allingment and
anisotropic viscosity are not included. The molecular field and the
passive stress are calculated by the functions

``` {.python language="Python"}
calc_molecular_field(self,Q)
calc_passive_stress_f(self,Q)
```

The linear and non-linear part of the evolution equation for $Q$,
eq. ([\[eq:EB\]](#eq:EB){reference-type="ref" reference="eq:EB"}) is
given as $$\begin{aligned}
    \omega(\nabla) &= \frac{K}{\gamma} \nabla^2 +\frac{AB}{\gamma}, \\
    N(Q) &= - \mathbf u\cdot \nabla Q + Q \Omega -\Omega Q - \frac{2A}{\gamma}Q^2_{kk}Q
\end{aligned}$$ The evolution of this is handled by the function

``` {.python language="Python"}
nem.evolve_nematic(self,number_of_steps,method='ETD2RK')
```

## Disipative dynamics

Note that if we set the velocity field to zero the dynamics become
$$\partial_t Q=  \frac{K}{\gamma} \nabla^2 Q_{ij} +\frac{A}{\gamma}(B - 2Q^2_{kk})Q_{ij}.$$
This is used to relax the initial system before starting the simulation.
The linear and nonlinear part of this equation are $$\begin{aligned}
    \omega(\nabla) = \frac{K}{\gamma} \nabla^2 +\frac{AB}{\gamma},  \\
    N(Q) = - \frac{2A}{\gamma}Q^2_{kk}Q.
\end{aligned}$$ An evolver for this dissipative dynamics is included as

``` {.python language="Python"}
nem.evolve_nematic_no_flow(self,number_of_steps,method='ETD2RK')
```

## The velocity field {#sec:nem_vel}

For a given orderparameter $Q$ the velocity field is calculated in
Fourier space using
eq.([\[eq:Stokes_nematic\]](#eq:Stokes_nematic){reference-type="ref"
reference="eq:Stokes_nematic"}). We start by finding an expression for
the pressure by taking the divergence of this eq.
([\[eq:Stokes_nematic\]](#eq:Stokes_nematic){reference-type="ref"
reference="eq:Stokes_nematic"}) and then using the incompressibility
condition giving $$\nabla^2 P = \nabla \cdot \mathbf F,$$ where
$\mathbf F =  \nabla \cdot \sigma^a(Q) +\nabla \cdot \sigma^p(Q)$ is the
active and passive forces. This is solved in Fourier space as
$$-k^2  {{P}_{\scriptscriptstyle  f}} =  i\mathbf k \cdot {{\mathbf F}_{\scriptscriptstyle  f}}.$$
The above equation can be inverted in order to find all the modes of the
pressure except the zero mode, i.e the pressure is determined up to a
constant. We set this constant to zero. Once the pressure is found we
obtain the velocity from
$$(\Gamma + \eta k^2){{\mathbf u}_{\scriptscriptstyle  f}} = - i\mathbf k {{P}_{\scriptscriptstyle  f}} + {{F}_{\scriptscriptstyle  f}}.$$
Note that when $\Gamma = 0$ we need to set the zero mode of the velocity
by hand. This is set to zero. The pressure and velocity are
calculated/updated by the two functions

``` {.python language="Python"}
calc_pressure_f(self) 
conf_u(self,Q)
```

Note that `calc_pressure_f` only returns the Fourier transform of the
pressure. The function `conf_u` updates both the velocity field `self.u`
and its Fourier transform `self.u_f`.

## Topological defects and active turbulence

Because of the head-tail symmetry if the nematic director the
topological defects in the nematic phase can have half integer winding
number. We can see this by maping the $Q$ tensor to a complex field.
This is done by writing the nematic director as
$\\vec{n} = \cos{\theta} \hat x + \sin{\theta} \hat y$, with
$\hat x /\hat y$ being the unit vectors in $x /y$ direction, and mapping
the $Q$ tensor, see
eq. ([\[eq:Q_tensor\]](#eq:Q_tensor){reference-type="ref"
reference="eq:Q_tensor"}), to the complex field
$$\psi = Q_{xx} +  iQ_{xy} = \frac{S}{2} e^{2 i\theta}.$$
Using the same arguments as for the BEC we find that the allowed winding
numbers $$k = \int_C \nabla \theta \cdot d\\vec{l} = 2\pi q$$ with
$q$ being a half-integer. The defects of lowest absolute charge is the
$\pm 1/2$ defects, which are depicted in fig.
[7.1](#fig:nem_dipole){reference-type="ref" reference="fig:nem_dipole"}.

![The nematic director (head-less vectors) around a defect dipole. The
$+1/2$ defect is marked with red, while the $-1/2$ defect is in blue.
](Figures/Dipole_nem.png){#fig:nem_dipole}

For tracing the defect nodes one can use the function

``` {.python language="Python"}
calc_vortex_nodes_nem(self, dt_Q=None,polarization = None)
```

If `dt_Q` is given this finds the defects velocity and if
`polarization ` is given the polarization of the $+1/2$ defects are
found. This polarization is given by $$\\vec{e}_+ = 
    \left( \frac{\nabla \cdot Q}{|\nabla \cdot Q|}\right)_{\\vec{r}= \\vec{r}_+}$$
where $\\vec{r}_+$ is the defects position. The field
$\\vec{e}_+$ can be found by the function

``` {.python language="Python"}
calc_defect_polarization_field(self)
```

Note that this function does not include the normalization to avoid
division by zero.


## Initial States

A homogeneous nematic with random noise can be implemented with the
function

``` {.python language="Python"}
conf_initial_condition_disordered(self, noise_strength=0.01)  
```

This gives a state where the nematogens are aligned with the $x$-axis.
The noise is in the angle $\theta$. If the activity is high enough this
state will after a while start to spontaneously generate topological
defects. In addition there is possible to insert defect dipoles using
the function

``` {.python language="Python"}
conf_insert_vortex_dipole(self, dipole_vector=None, dipole_position=None) 
```

which works similarly as the one implemented for the BEC. This function
can be used either to initialize a homogeneous state with a dipole, or
it can be used to insert a dipole into an already existing nematic.

## Spatially varying activity

The activity $\alpha$ is can be spatially varying. This can be used to
make active channels as shown in
fig. [7.2](#fig:active_channel){reference-type="ref"
reference="fig:active_channel"}

![Illustration of an active channel. $\text{width}= 20$ and
$d = 2$.](Figures/FiguresOfSetup/spatial_alpha.png){#fig:active_channel}

This simple channel with the activity $\alpha_0$ inside and $\alpha = 0$
outside is included as the function

``` {.python language="Python"}
conf_active_channel(self,width,d=7)
```

Which sets the activity to
$$\alpha = \alpha_0 \left[1 -1/2 \left(\tanh([x-w/2]/d) -\tanh([x+w/2]/d) \right)\right].$$
$w$ is here the width of the channel, $\alpha_0$ is the activity in the
channel and $d$ is the width of the interface between the channel and
the bulk. More complicated structures can be created if one wish.


## Topological defects

Topological defects in nematic liquid crystals are called disclinations and are characterized the orientation of the rod-like particles having rotated after following a path around the dislocation. 
From [schimmingKinematicsDynamicsDisclination2023](References.md),

$$
D_{\gamma i} = \epsilon_{\gamma \mu \nu} \epsilon_{ikl} \partial_k Q_{\mu \alpha} \partial_l Q_{\nu \alpha}.
$$

In two dimensions, where $\mathbf n = (\cos \theta,\sin \theta)$, we have 

$$
Q = 
S\begin{pmatrix}
\cos \theta \cos \theta - \frac{1}{2} & \cos \theta \sin \theta \\
\cos \theta \sin \theta & \sin \theta \sin \theta - \frac{1}{2}\\
\end{pmatrix}
=
\frac{S}{2}
\begin{pmatrix}
\cos (2\theta) &  \sin(2\theta) \\
\sin(2\theta) & - \cos (2\theta)\\
\end{pmatrix},
$$

$$
= 
\frac{1}{2}
\begin{pmatrix}
\psi_1 &  \psi_2 \\
\psi_2 & - \psi_1\\
\end{pmatrix},
$$

where $\psi_1,\psi_2$ are the components of an $\mathcal D^2$ order parameter.
We get only one component of $D_{\gamma i}$, which is $D_{33}$ which is

$$
D_{33} = \epsilon_{\mu \nu} \epsilon_{kl} \partial_k Q_{\mu \alpha} \partial_l Q_{\nu \alpha}
$$
$$
=
\epsilon_{\mu \nu} \epsilon_{kl} \partial_k Q_{\mu 1} \partial_l Q_{\nu 1}
+
\epsilon_{\mu \nu} \epsilon_{kl} \partial_k Q_{\mu 2} \partial_l Q_{\nu 2}
$$

We have $Q_{\mu1} = \frac{1}{2} \psi_\mu$ and $Q_{\mu 2} = \frac{1}{2} \epsilon_{\mu q} \psi_q$, so 

$$
D_{33}
=
\frac{1}{4} \epsilon_{\mu \nu} \epsilon_{kl} (\partial_k \psi_\mu )(\partial_l \psi_\nu)
+
\frac{1}{4} \epsilon_{\mu \nu} \epsilon_{kl} (\partial_k  \epsilon_{\mu q} \psi_q) (\partial_l \epsilon_{\nu r} \psi_r)
$$
And using that 

$$\epsilon_{\mu \nu} \epsilon_{\mu q} \epsilon_{\nu r} = \epsilon_{qr},
$$

we get

$$
D_{33} = \frac{1}{2} \epsilon_{\mu \nu} \epsilon_{kl} (\partial_k \psi_\mu )(\partial_l \psi_\nu).
$$

This is the same determinant as we would get using the coarse grain density of [skogvollUnifiedFieldTheory2023](References.md), only with $\psi_0$, so, the disclination density should be 

$$
\rho_{\gamma i} = \frac{1}{\pi S_0^2} D_{\gamma i} 
$$

In three dimenstions, the story is more complicated, because we have a tensor $\rho_{\gamma i}$. 
This tensor contains two pieces of information, namely which direction the disclination is pointing, and around which axis $\boldsymbol \Omega$, near the disclination, the rods are rotating. 
In that way, it is similar to a dislocation density in a crystal structure, only that it allows for the orientation the "Burgers vector" to be any direction. 
It can probably be written like this

$$
\rho_{\gamma i} = \Omega_\gamma t_{i} ,
$$
where $\boldsymbol t$ is a unit vector. 
From this, we see that 

$$
|\rho|^2 = \rho_{\gamma i} \rho_{\gamma i} = |\boldsymbol \Omega|^2,
$$

so $\sqrt{|\rho|^2}$ is the quantity we should integrate to find the nodes of the defects.

From [schimmingKinematicsDynamicsDisclination2023](References.md), we have 

$$
t_i \delta^{(2)}(\mathbf r_{\perp}) = \delta^{(2)}(\mathbf Q_\perp) \Omega_\gamma D_{\gamma i}
$$

replacing the delta function, which we may generalize to 

$$
\mathbf \rho = \Omega_{\gamma} \frac{1}{\pi Q_{\perp}^2} D_{\gamma i}.
$$

