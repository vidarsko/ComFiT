# Class: Nematic Liquid Crystal

A liquid crystal is a state of matter between a solid and a liquid.

In this class, we simulate an active nematic liquid crystal using [framework].

```python
file: comfit/models/nematic_liquid_crystal.py 
class: NematicLiquidCrystal
```

## Variables and parameters

The primary variables are the symmetric traceless tensor $Q$ and the
velocity field $\\mathbf{u}$

```python
nem.Q 
nem.u
```

The $Q$ tensor is given by the nematic director as

$$
Q_{ij} = S (n_i n_j - \frac{1}{d} \delta_{ij})
$$

in $d$ dimensions.
To take advantage of its
symmetric nature we have saved $Q$ as a vector field, which in two and three
dimensions takes the forms

$$
\begin{aligned}
\mathbf{Q} &= [ Q_{xx},Q_{xy} ] \\
\mathbf{Q} &= [ Q_{xx},Q_{xy},Q_{xz},Q_{yy},Q_{yz}]
\end{aligned}
$$

respectivly.
We can translate from tensor indexes to the right value stored in the vector by
using the function

```python
get_sym(self,Q,i,j)
```

This returns the element $Q_{ij}$ of a symmetric and traceless tensor
field. In addition to this we also have the function

```python
get_anti_sym(self,omega,i,j)
```

so that we can optimally store the antisymetric tensors as well. In two
dimensions these only have one independent component, which is stored as
a scalar field,
while in three dimensions it is stored as

$$
\mathbf \Omega = [\Omega_{xy}, \Omega_{xz}, \Omega_{yz}]
$$

In order to calculate the director field $\mathbf n$ and the amount of order $S$ in two dimensions we use that we can map the orderparameter to the complex field $\psi = Q_{xx} +  iQ_{xy} =Se^{2i\theta}/2$, where $\theta$  is the angle of the director field.
In three dimensions we use that $S$ is given by the largest eigenvalue as $S = 3\lambda/2$ with the director being the coresponding eigenvector [^Schimming2022Thesis].
This is taken care of in the function

```python
calc_order_and_director(self)
```

Wich returns `S, n`, where `n` is the director field.
Note that a nematic liquid crystall can be biaxial and given as

$$
Q_{ij} = S (n_i n_j - \frac{1}{3} \delta_{ij}) + P (m_i m_j -l_i l_j)
$$

where $P$ is given by the difference between the smallest eigenvalues and $\mathbf m$ and $\mathbf l$ is the corresponding eigenvectors.

## Model

We model the active nematic using a set of coupled differential
equations, namely the Edvard-Beris equation coupled to the Stokes
equation
[^marchetti2013hydrodynamics] [^genkin2017topological] [^nejad2020memory] [^angheluta2021role]

$$
\begin{aligned}
\partial_t Q + \mathbf u\cdot \nabla Q +\Omega Q -Q \Omega &=\gamma^{-1}H,
\\
(\Gamma- \eta \nabla^2 )\mathbf u &= -\nabla P + \nabla \cdot \sigma^a(Q) + \nabla \cdot \sigma^p, \\
    \nabla \cdot \mathbf u &= 0.
\end{aligned}
$$

Here $2\Omega_{ij} = \partial_i u_j - \partial_j u_i$ is
the vorticity tensor, $P$ is the pressure, $\gamma$ is the rotational friction coefficient, $\sigma^p$ is the passive stress, $\Gamma$ is friction with a substrate, $\eta$ is viscosity and the active stress is given by $\sigma^a = \alpha Q$.
The vorticity tensor is calculated by the
function

```python
calc_vorticity_tensor(self)
```

Note that the velocity has to be updated before this function is called.
The calculation of the pressure and velocity is described furhter down.
Since the active stress is simply proportional to $Q$ we have not included any
function to calculate it, but calculate the force directly with the function

```python
calc_active_force_f(self,Q)
```

The molecular field $H$ is given as

$$H_{ij} =  -\frac{\delta \mathcal{F}}{\delta Q_{ij}} + \frac{\delta_{ij}}{d} \text{Tr}\left(\frac{\delta F}{\delta Q}\right)
$$

The last term is here to make it trace less.
For the free energy we use

$$\mathcal F = \int \left( K |\nabla Q|^2 - \frac{A}{2} \left[ B \text{Tr}(Q^2) -\text{Tr}(Q^2)^2   \right] -\frac{C}{3}\text{Tr}(Q^3) \right),
$$

where it is assumed that there is a single Frank elastic constant $K$.
Note that the last term only exists in three dimensions since the trace of $Q^3$ vanishes in two dimensions.
The molecular field is then given as

$$
 H_{ij} =  K \nabla^2 Q_{ij} + A(B - 2Q^2_{kk})Q_{ij} +
 \begin{cases}
  0, & \text{dim} = 2 \\
 C Q^2_{ij} - \frac{C}{3}Q^2_{kk} \delta_{ij}, & \text{dim} = 3
\end{cases}
$$

For the passive stress we use
$$\sigma_{ij} = Q_{ik}H_{kj}-H_{ik}Q_{kj} - K (\partial_i Q_{kl})(\partial_jQ_{kl}).$$
Note that we use the convention that the dot product between a vector
and a tensor contracts the last component when calculating the
divergence of this. The first two terms is the asymmetric stress, while
the second term is the Ericksen stress. Terms due to flow allingment and
anisotropic viscosity are not included. The molecular field and the
passive stress  are calculated by the functions

```python
calc_molecular_field(self,Q)
calc_passive_stress_f(self,Q)
```

The linear and non-linear part of the evolution equation for $Q$,
eq. is
given as

$$
\begin{aligned}
\omega(\nabla) &= \frac{K}{\gamma} \nabla^2 +\frac{AB}{\gamma}, \\
N(Q,t) &= - \mathbf u\cdot \nabla Q + Q \Omega -\Omega Q - \frac{2A}{\gamma}Q^2_{kk}Q +\begin{cases}
  0, & \text{dim} = 2 \\
 C Q^2 - \frac{C}{3}Q^2_{kk} I, & \text{dim} = 3
\end{cases}
\end{aligned}
$$

The evolution of this is handled by the function

```python
nem.evolve_nematic(self,number_of_steps,method='ETD2RK')
```

## Disipative dynamics

Note that if we set the velocity field to zero the dynamics become
$$\partial_t Q=  \frac{K}{\gamma} \nabla^2 Q_{ij} +\frac{A}{\gamma}(B - 2Q^2_{kk})Q_{ij}.$$
This is used to relax the initial system before starting the simulation.
The linear and nonlinear part of this equation are

$$
\begin{aligned}
    \omega(\nabla) = \frac{K}{\gamma} \nabla^2 +\frac{AB}{\gamma},  \\
    N(Q) = - \frac{2A}{\gamma}Q^2_{kk}Q +\begin{cases}
  0, & \text{dim} = 2 \\
 C Q^2 - \frac{C}{3}Q^2_{kk} I, & \text{dim} = 3
\end{cases}.
\end{aligned}
$$

An evolver for this dissipative dynamics is included as

```python
nem.evolve_nematic_no_flow(self,number_of_steps,method='ETD2RK')
```

## The velocity field {#sec:nem_vel}

For a given orderparameter $Q$ the velocity field is calculated in
Fourier space using
eq (??).
We start by finding an expression for
the pressure by taking the divergence of this eq. (??) and then using the incompressibility
condition giving $$\nabla^2 P = \nabla \cdot \mathbf F,$$ where
$\mathbf F =  \nabla \cdot \sigma^a(Q) +\nabla \cdot \sigma^p(Q)$ is the
active and passive forces. This is solved in Fourier space as

$$
-k^2  {{P}_{\scriptscriptstyle  f}} =  i\mathbf k \cdot {{\mathbf F}_{\scriptscriptstyle  f}}.
$$

The above equation can be inverted in order to find all the modes of the
pressure except the zero mode, i.e the pressure is determined up to a
constant. We set this constant to zero. Once the pressure is found we
obtain the velocity from

```math
(\Gamma + \eta k^2){{\mathbf u}_{\mathfrak  f}} = - i\mathbf k {{P}_{\mathfrak  f}} + {{F}_{\mathfrak  f}}.
```

Note that when $\Gamma = 0$ we need to set the zero mode of the velocity
by hand. This is set to zero. The pressure and velocity are
calculated/updated by the two functions

```python
calc_pressure_f(self) 
conf_u(self,Q)
```

Note that `calc_pressure_f` only returns the Fourier transform of the
pressure. The function `conf_u` updates both the velocity field `self.u`
and its Fourier transform `self.u_f`.

## Minimum of the free energy

When starting a simulation it is often interesting to start from a configuration that is the minimum of the free energy pluss some perturbations or with a vortex dipole/fillament.
From the free energy we see that the minimum energy is given by a homogeneous nematic, and it is inedependent of the direction the nematogens are pointing.
Assuming that the unitvector $\mathbf n$ is homogeneous we can rewrite the free energy in terms of the parameter $S$.

### In two dimmensions

the free energy is only given by powers of $\text{Tr}(Q^2)$ which is $S^2/2$ in terms of $S$.
The free-energy is therfore for a homogeneous two dimentional nematic given as

$$
\mathcal F =  \int \left( - \frac{A}{2} \left[ B \frac{S^2}{2} -\frac{S^4}{4}   \right] \right).
$$

The minimum of this is given as $S =\sqrt{B}$ when $B >0$ and $S = 0$ if $B<0$.

### In three dimensions

In three dimensions we have that $\text{Tr}(Q^2) = 2 S^2/3$ and $\text{Tr}(Q^3)= 2 S^3/9$.
Using this we find that there are a local minima at

$$
S = \frac{1}{8}\frac{C}{A} + \frac{1}{2} \sqrt{\frac{C^2}{16 A^2} + 3 B}
$$

when $B > -3 C^2/(16A^2)$.

## Topological defects and active turbulence

Because of the head-tail symmetry if the nematic director the
topological defects in the nematic phase can have half integer winding
number. We can see this by maping the $Q$ tensor to a complex field.
This is done by writing the nematic director as
$\\mathbf{n} = \cos{\theta} \hat x + \sin{\theta} \hat y$, with
$\hat x /\hat y$ being the unit vectors in $x /y$ direction, and mapping
the $Q$ tensor, see
eq. (??), to the complex field

$$
\psi = Q_{xx} +  iQ_{xy} = \frac{S}{2} e^{2 i\theta}.
$$

Using the same arguments as for the BEC we find that the allowed winding
numbers $$k = \int_C \nabla \theta \cdot d\\mathbf{l} = 2\pi q$$ with
$q$ being a half-integer. The defects of lowest absolute charge is the
$\pm 1/2$ defects, which are depicted below.

![Liquid crystal disclination dipole](images/nematic_liquid_crystal_disclination_dipole.png#only-light)
![Liquid crystal disclination dipole](images/nematic_liquid_crystal_disclination_dipole-colorinverted.png#only-dark)

*Liquid crystal disclination dipole:* The nematic director (head-less vectors) around a defect dipole. The $+1/2$ defect is marked with red, while the $-1/2$ defect is in blue.

For tracing the defect nodes one can use the function

```python
calc_vortex_nodes_nem(self, dt_Q=None,polarization = None)
```

If `dt_Q` is given this finds the defects velocity and if `polarization` is given the polarization of the $+1/2$ defects are
found.
This polarization is given by

$$
\mathbf{e}_+ = \left( \frac{\nabla \cdot Q}{|\nabla \cdot Q|}\right)_{\mathbf{r}= \mathbf{r}_+}
$$

where $\mathbf{r}_+$ is the defects position.
The field $\mathbf{e}_+$ can be found by the function

```python
calc_defect_polarization_field(self)
```

Note that this function does not include the normalization to avoid
division by zero.

## Initial States

A homogeneous nematic with random noise can be implemented with the
function

```python
conf_initial_condition_disordered(self, noise_strength=0.01)  
```

This gives a state where the nematogens are aligned with the $x$-axis.
The noise is in the angle $\theta$. If the activity is high enough this
state will after a while start to spontaneously generate topological
defects. In addition there is possible to insert defect dipoles using
the function

```python
conf_insert_vortex_dipole(self, dipole_vector=None, dipole_position=None) 
```

which works similarly as the one implemented for the BEC. This function
can be used either to initialize a homogeneous state with a dipole, or
it can be used to insert a dipole into an already existing nematic.

## Spatially varying activity

The activity $\alpha$ is can be spatially varying. This can be used to
make active channels as shown in the following figure.

![Nematic liquid crystal active channel](images/nematic_liquid_crystal_active_channel.png#only-light)
![Nematic liquid crystal active channel](images/nematic_liquid_crystal_active_channel-colorinverted.png#only-dark)

*Nematic liquid crystal active channel:* Illustration of an active channel. `width= 20` and $d = 2$.

This simple channel with the activity $\alpha_0$ inside and $\alpha = 0$
outside is included as the function

```python
conf_active_channel(self,width,d=7)
```

Which sets the activity to
$$\alpha = \alpha_0 \left[1 -1/2 \left(\tanh([x-w/2]/d) -\tanh([x+w/2]/d) \right)\right].$$
$w$ is here the width of the channel, $\alpha_0$ is the activity in the
channel and $d$ is the width of the interface between the channel and
the bulk. More complicated structures can be created if one wish.

## Three dimensions

In three dimensions, the $Q$-tensor can be characterized as

$$
Q_{ij} = S (n_i n_j - \frac{1}{3} \delta_{ij} ) + P (m_i m_j - l_i l_j),
$$

where $(\mathbf n, \mathbf m, \mathbf l )$ is an orthornmal triad.
It is five parameters: $S$, the two defining angles of $\mathbf n$, $P$ and the angle of $\mathbf m$ in the plane orthogonal to $\mathbf n$.
$S$ can always be determined from the highest eigenvalue $\lambda_{\textrm{max}}$ of $Q$ by [^Schimming2022Thesis]

$$
S = \frac{3}{2} \lambda_{\textrm{max}}
$$

## Topological defects

Topological defects in nematic liquid crystals are called disclinations and are characterized the orientation of the rod-like particles having rotated after following a path around the dislocation.
From [^schimming2023kinematics],

$$
D_{\gamma i} = \epsilon_{\gamma \mu \nu} \epsilon_{ikl} \partial_k Q_{\mu \alpha} \partial_l Q_{\nu \alpha}.
$$

In two dimensions, where $\mathbf n = (\cos \theta,\sin \theta)$, we have

$$
Q = S\begin{pmatrix}
\cos \theta \cos \theta - \frac{1}{2} & \cos \theta \sin \theta \\
\cos \theta \sin \theta & \sin \theta \sin \theta - \frac{1}{2}\\
\end{pmatrix}
= \frac{S}{2}
\begin{pmatrix}
\cos (2\theta) &  \sin(2\theta) \\
\sin(2\theta) & - \cos (2\theta)\\
\end{pmatrix},
$$

$$
= \frac{1}{2}
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

```math
= \epsilon_{\mu \nu} \epsilon_{kl} \partial_k Q_{\mu 1} \partial_l Q_{\nu 1}
+ \epsilon_{\mu \nu} \epsilon_{kl} \partial_k Q_{\mu 2} \partial_l Q_{\nu 2}
```

We have $Q_{\mu1} = \frac{1}{2} \psi_\mu$ and $Q_{\mu 2} = \frac{1}{2} \epsilon_{\mu q} \psi_q$, so

```math
D_{33} = \frac{1}{4} \epsilon_{\mu \nu} \epsilon_{kl} (\partial_k \psi_\mu )(\partial_l \psi_\nu)
+ \frac{1}{4} \epsilon_{\mu \nu} \epsilon_{kl} (\partial_k  \epsilon_{\mu q} \psi_q) (\partial_l \epsilon_{\nu r} \psi_r)
```

And using that

$$
\epsilon_{\mu \nu} \epsilon_{\mu q} \epsilon_{\nu r} = \epsilon_{qr},
$$

we get

$$
D_{33} = \frac{1}{2} \epsilon_{\mu \nu} \epsilon_{kl} (\partial_k \psi_\mu )(\partial_l \psi_\nu).
$$

This is the same determinant as we would get using the coarse grain density of [^skogvoll2023Topological], only with $\psi_0$, so, the disclination density should be

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

From Ref.[^schimming2023kinematics], we have

$$
t_i \delta^{(2)}(\mathbf r_{\perp}) = \delta^{(2)}(\mathbf Q_\perp) \Omega_\gamma D_{\gamma i}
$$

replacing the delta function, which we may generalize to

$$
\mathbf \rho = \Omega_{\gamma} \frac{1}{\pi (S_0-P_0)^2} D_{\gamma i}.
$$

### Inserting topological defects

How do we insert topological defects of a given character?

We can generate an initial state of $Q_{ij}$ by writing

$$
Q_{ij} = S_0 \left (\frac{1}{2} n_i n_j - \frac{1}{d} \delta_{ij} \right ),
$$

and then simply impose an orientation field corresponding to an angle field on the $\mathbf n$ fields.

[^Schimming2022Thesis]:Schimming, C. D. (2022). Theoretical and Computational Methods for Mesoscopic Textures in Nematic Liquid Crystals with Anisotropic Elasticity. PhD Thesis. The University of Minnesota. [https://hdl.handle.net/11299/241713](https://hdl.handle.net/11299/241713)
[^marchetti2013hydrodynamics]: Marchetti, M. C., Joanny, J-F., Ramaswamy, S., Liverpool, T. B., Prost, J., and Rao, M. and Simha, R. A. (2013). Hydrodynamics of soft active matter. Reviews of Modern Physics. 85, 3, 1143. [https://doi.org/10.1103/RevModPhys.85.1143](https://doi.org/10.1103/RevModPhys.85.1143)
[^genkin2017topological]: Genkin, M. M., Sokolov, A., Lavrentovich, O. D. and Aranson, I. S. (2017). Topological defects in a living nematic ensnare swimming bacteria. Physical Review X. 7, 1,011029. [https://doi.org/10.1103/PhysRevX.7.011029](https://doi.org/10.1103/PhysRevX.7.011029)
[^nejad2020memory]: Nejad, M. R., Doostmohammadi, A. and Yeomans, J. M. (2021). Memory effects, arches and polar defect ordering at the cross-over from wet to dry active nematics. Soft Matter. 17, 9, 2500-2511. [https://doi.org/10.1039/D0SM01794A](https://doi.org/10.1039/D0SM01794A)
[^angheluta2021role]:Angheluta, L., Chen, Z., Marchetti, M. C. and Bowick, Mark J. (2021). The role of fluid flow in the dynamics of active nematic defects. New Journal of Physics. 23, 3, 033009. [https://doi.org/10.1088/1367-2630/abe8a8](https://doi.org/10.1088/1367-2630/abe8a8)
[^schimming2023kinematics]: Schimming, C. D. and Viñals, J. (2023). Kinematics and dynamics of disclination lines in three-dimensional nematics. Proceedings of the Royal Society A. 479, 2273, 20230042. [https://doi.org/10.1098/rspa.2023.0042](https://doi.org/10.1098/rspa.2023.0042)
[^skogvoll2023Topological]: Skogvoll, V., Rønning, J., Salvalaglio, M., Angheluta, L. (2023). A unified field theory of topological defects and non-linear local excitations. npj Comput Mater, 9, 122. [https://doi.org/10.1038/s41524-023-01077-6](https://doi.org/10.1038/s41524-023-01077-6)
