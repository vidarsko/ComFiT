# Class: Bose-Einstein Condensate

There are two types of particles in the world: fermions and bosons.
Whereas fermions can never occupy the same quantum state due to the Pauli Exclusion principle, the same is not true for bosons.
A Bose-Einstein condensate (BEC) is a state of matter consisting of ultra-cold bosons that undergo a phase transition at a low critical temperature in which most bosons occupy the ground state of the system.
It was theorized by Einstein and Bose in the 1920s as a new state of matter and produced for the first time in 1995 by Eric Cornell and Carl Wieman[^andersonObservationBoseEinsteinCondensation1995].

In this class, we simulate a Bose-Einstein condensate in 1, 2 and 3 dimensions using the Gross-Pitaevski equation (GPE).

```python
file: comfit/models/bose_einstein_condensate.py 
class: BoseEinsteinCondensate
```

## Variables and parameters

The primary field is the complex wave function $\psi$

```python
bec.psi 
```
The Bose-Einstein condensate takes the same keyword arguments as the BaseSystem in addition to

| Keyword | Definition | Default value|
|---------|------------|--------------|
| `gamma`  | dissipative factor | $0.01$ |

## Model

The BEC is in the mean field regime described by the GPE[^dalfovo1999theory] [^kevrekidis2007emergent]. 
This is a non-linear Schrödinger equation which reads

$$
i\hbar \partial_t\psi = \left[-\frac{\hbar^2}{2m} \nabla^2+ V_{ext} -\mu +g|\psi|^2 \right]\psi.
$$

Here, $\mu$ is the chemical potential, $m$ is the mass of the bosons, $g$ is an interaction parameter and $\hbar$ is the Planc constant. 
$\psi$ is the wave function describing the condensate phase and $V_{ext}$ is an external potential. 
The GPE can be obtained from the variational principle [^kevrekidis2007emergent][^pitaevskiiBook]

$$
\mathfrak i \hbar \partial_t \psi = \frac{\delta K}{\delta \psi^*}
$$

with the Hamiltonian

$$
K = \int d \mathbf r \left[\frac{\hbar^2}{2m}|\nabla\psi|^2 +(V_{ext} -\mu)|\psi|^2 +\frac{g}{2}|\psi|^4 \right].
$$

We introduce dimensionless units for length $\xi = \hbar/\sqrt{m\mu}$, time $\tau = \xi/c$ and energy $E=\eta$ and rescaling the wave function to $\psi \rightarrow \sqrt{\frac{g}{\mu}}\psi$, in addition we include a dissipative factor $\gamma$.
This results in the dGPE on dimensionless form as
[^gardiner2003stochastic] [^rooney2012stochastic] [^bradley2012energy] [^skaugenUnifiedPerspectiveTwodimensional2018]

<!-- $$
\mathfrak i \partial_t \psi = (1-\mathfrak i\gamma) \left[-\frac{1}{2}\nabla^2 + V_{ext} -1 +|\psi|^2 \right]\psi.
$$

$$
\partial_t \psi =-\mathfrak i (1-\mathfrak i\gamma) \left[-\frac{1}{2}\nabla^2 + V_{ext} -1 +|\psi|^2 \right]\psi.
$$

$$
\partial_t \psi =-(\mathfrak i+\gamma) \left[-\frac{1}{2}\nabla^2 + V_{ext} -1 +|\psi|^2 \right]\psi.
$$ -->

!!! equation "The danmped Gross-Pitaevski equation (`evolve_dGPE`)"
    $$
    \partial_t \psi =-(\mathfrak i+\gamma) \left[-\frac{1}{2}\nabla^2 + V_{ext} -1 +|\psi|^2 \right]\psi.
    $$

So we see that we can split into its linear and non-linear parts as

$$
\omega = (\mathfrak{i}+\gamma) (1 + \frac{1}{2}\nabla^2), \quad N = - (\mathfrak{i} + \gamma) (V_{ext} + |\psi|^2)\psi
$$

The Hamiltonian and its density in dimensionless units are calculated by the functions

```python
bec.calc_hamiltonian(self)
bec.calc_hamiltonian_density(self)
```

## Approximation of Ground States

In simulations, it is often convenient to start in a configuration that is close to the ground state. 
We can estimate this ground state by noticing that the GPE dissipates energy when it is evolved in imaginary time $t \rightarrow it$
[^minguzzi2004numerical] [^kevrekidis2007emergent].
We refer to this as *evolving the dGPE in imaginary time*.
Given an external potential $V_{ext}$ we can therefore find an approximation to the ground state by starting with a guess and then removing energy from the guessed state by evolving the equations in imaginary time.

To get a guess of the ground state for a given potential $V_{ext}$ we use the Thomas-Fermi approximation [^dalfovo1999theory] [^kevrekidis2007emergent], where we assume that $\psi$ is slowly varying so that we can neglect the Laplacian term.
Looking for stationary solutions to the dGPE we obtain the equation

$$
0 = (V_{ext} - 1 +|\psi|^2 )\psi.
$$

This has two solutions, $\psi = 0$
and

$$
|\psi|^2 = 1-V_{ext}.
$$

In the case where $V_{ext} > 1$ there is only one possibility, namely $\psi = 0$. 
In the case of $V_{ext}  < 1$, both the stationary solutions exist, so we need to evaluate their energy. 
We can do this by considering the Hamiltonian

$$
K_{TF} = \int d\mathbf r \left[(V_{ext}  -1)|\psi|^2 + \frac{1}{2} |\psi|^4 \right].
$$

Inserting the anzats $\psi =0$ we see that this vanishes, while if we insert the anzats $|\psi|^2 = 1-V_{ext}$, we get something negative. 
We therefore conclude that the Thomas-Fermi ground state is given as

$$
\psi = \begin{cases}
     0 & \text{if} \quad V_{ext}  > 1 \\
     \sqrt{1 -V_{ext}}  & \text{if}\quad V_{ext}  < 1
    \end{cases}
$$

This ground state can be initialized as

```python
bec.conf_initial_condition_Thomas_Fermi(self)
```

Once this guess is initialized we can propagate the wave function in imaginary time using the function

```python
bec.evolve_relax_BEC(self,number_of_steps,method='ETD2RK')
```

Note that the potential has to be a constant during this evolution.

## Potentials

A lot of the interesting dynamics of a BEC come about when it is interacting with an external potential. 
This is included as the function

```python
bec.V_ext(t) 
```

By default, this is assumed to be time-independent and returns the field

```python
bec.V_0 
```

The potential can be changed by the function

```{.python language="Python"}
bec.conf_external_potential(self, V_ext, additive=False)
```

which can be used both to set it as a function or to set it as a constant potential depending on whether `V_ext` is a function, a constant or an numpy array. 
If `additive =True` one add the constant `V_ext` to the allready existing potential.
If `V_ext` is a function it needs to be on the form

```{.python language="Python"}
def V(t)
     ...
     return ...
```

The evolver will then evaluate it using the `bec.time` variable which is updated on the run.
An example using a Gaussian stirring potential is provided in the example folder.

To make life easier we have provided a couple of popular potentials.
The harmonic potential is provided in the function

```python
bec.set_harmonic_potential(self,R_tf)
```

Here $R_{tf}$ is the Thomas-Fermi radius [^kevrekidis2007emergent], and the harmonic potential takes the form $$V_H = \frac{r^2}{R_{tf}^2}$$.
The Gaussian potential is provided through the function

```python
bec.gaussian_stirring_potential(self,size,strength,position)
```

giving the potential

$$
V_g = g e^{-|\mathbf{r} - \mathbf{r}_p|^2/\sigma^2}.
$$

Here
$g$ is the strength, $\sigma$ is the size and $\mathbf{r}_p$ is the position.

Much of the interesting physics happens when the potential is time-dependent. 
For this one can define a function as

```python
def Func():
    ...
```

and update the external potential to this function by calling

## Hydrodynamics

The dGPE can be transformed to a hydrodynamic description of the BEC.
The first step in doing this is to introduce the Madelung transformation $\psi = \sqrt{\rho} e^{i\theta}$, where $\rho = |\psi|^2$ is the superfluid density [^kevrekidis2007emergent]. 
For $\gamma = 0$ this density is conserved and satisfies the conservation equation

$$
\partial_t \rho + \nabla\cdot \mathbf{J}_s = 0,
$$

with the superfluid current given by

$$
\mathbf{J}_s = \Im(\psi^* \nabla \psi) = \rho \nabla \theta = \rho \mathbf{v}_s.
$$

Here the superfluid velocity is introduced as
$\mathbf{v}_s = \nabla \theta$. 
The superfluid current can be found with the function

```python
bec.calc_superfluid_current(self)   
```

We can also put the Madelung transformation into the Hamiltonian to get [^bradley2012energy] [^nore1997kolmogorov]

$$
K = \int d \mathbf r \left[\frac{1}{2}\rho v_s^2 +\frac{1}{8} \frac{|\nabla \rho|^2}{\rho} + (V_{ext}-1)\rho +\frac{1}{2}\rho^4 \right].
$$

The first term here is the kinetic energy of the condensate.
To calculate this it is convenient to introduce the density weighted velocity $\mathbf{u} = \sqrt{\rho}\mathbf{v}_s$[^bradley2012energy].
This have the advantage of not being singular at the centre of the topological defects. 
Using this we can write the kinetic energy as $$E_k =\int d \mathbf r \frac{1}{2}u^2.$$ 
The density weighted velocity and the kinetic energy can be calculated by the functions

```python
bec.calc_velocity(self)
bec.calc_kinetic_energy(self)
```

Further, if we insert the Madelung transformation into the dGPE and do some work we can map it into the Navier-Stokes equations [^kevrekidis2007emergent] [^bradley2012energy]

$$
\begin{aligned}
    \partial_t \rho + \nabla\cdot(\rho \mathbf v) =2\gamma \rho (1-V_{eff}),
    \\
    \partial_t \mathbf v-\mathbf v\cdot\nabla \mathbf v =  - \nabla(V_{ext} +\rho) +\frac{\gamma}{2}\nabla^2 \mathbf v.
\end{aligned}
$$

Notice that the condensate density is only conserved when $\gamma = 0$.

## Forces on external potential

To study how the BEC interacts with impurities, one can model the impurity as a Gaussian potential and measure the forces that are acting
on it [^ronning2020classical] [^astrakharchik2004motion] [^pinsker2017gaussian].
From the Ehrenfest theorem the forces on the condensate from the stirrer is given as

$$
\mathbf{F} = -\langle \nabla V_{ext}\rangle,
$$

which means that the force acting on the stirrer is $\mathbf{F} = \langle \nabla V_{ext}\rangle$. 
Written explicitly, this is

$$
\mathbf{F} = \int d\mathbf r |\psi^2| \nabla V_{ext} = -\int d\mathbf r V_{ext}\nabla|\psi^2|,
$$

and is calculated by the function

```python
bec.calc_force_on_external_potential(self)
```

Note that this calculates the average force on all of the external potential.

## Quantized vortices

The topological defects in the BEC take the form of quantized vortices.
This is because the velocity is given by the gradient of the phase $\mathbf{v}_s = \nabla \theta$. 
This has the consequence that the circulation $$\Gamma = \oint d\mathbf{l} \cdot \mathbf{v} = \oint d\mathbf{l} \cdot \nabla \theta = \int d\theta$$
is quantized in units of $2\pi$, since the field $\psi$ is single-valued.
In two dimensions this vortex is a point defect, while in three dimensions it is a line defect.
For the BEC this method is given by the function

```python
bec.calc_vortex_nodes( dt_psi=None)
```

which finds the vortex nodes in two or three dimensions. 
If $\partial_t \psi$ is put in as the field `dt_psi`, then the velocity of the defects are also found.

The defects can be created in multiple ways. In addition to putting them in by hand, one can create them by stirring the condensate with a potential or by relaxing a disordered initial state.

## Comoving frame

When studying a BEC that is stirred by a potential it is in some cases convenient to consider the potential as stationary with the BEC flowing around. 
This can be done by boosting the dGPE so that it reads

$$
\partial_t \psi = \mathbf V_p \cdot \nabla \psi +(\mathfrak i+\gamma) (1+\frac{1}{2}\nabla^2) \psi - (\mathfrak i + \gamma) (\mathcal U + |\psi|^2)\psi,
$$

where $\mathbf{V_p}$ is the velocity of the boost. Note that a Galilean boost of the GPE is often accompanied by a phase shift of the wave
function
$\psi \rightarrow \psi \exp{(\mathfrak i\mathbf V_p \cdot \mathbf r + \frac i 2 V_p^2 t)}$, which transforms the superfluid velocity to the new reference frame[^Pismen], leaving the GPE unchanged after the Galilean transformation.
However, the equation with $\gamma \neq 0$ is not Galilean invariant and we therefore do not include the phase factor in the transformation.
This has the consequence that the superfluid velocity obtained from $\psi$ is measured in the lab frame and not the comoving frame[^ronning2020classical].

To reduce the recycling of excitations into the incoming flow, we introduce a buffer region around the computational domain where $\gamma$ is large, similar to the one used in Ref. [^reeves2015identifying]. 
The dissipative factor becomes a function of space and is given by

$$
\gamma(\mathbf{r}) = \max[\gamma_x(x),\gamma_y(y),\gamma_z(z)]
$$ 

in three dimensions and

$$
\gamma(\mathbf{r}) =
\max[\gamma_x(x),\gamma_y(y)
$$ 

in two dimensions. Here

$$
\begin{aligned}
\gamma_x(x)= \frac{1}{2}\big(2 + \tanh{[(x-w_x)/d]}
-\tanh{[(x+w_x)/d]}\big) + \gamma_0.
\end{aligned}
$$

and similar for $\gamma_y(y)$ and $\gamma_z(z)$.
The constant $\gamma_0$ is the value in the bulk, $d$ is the size of the interface between the bulk and the buffer and $w_x$ is the distance from the centre of the domain to the buffer in the $x$ direction.
Note that when $\gamma$ is a function of space we can no longer put it into the non-linear differential operator $\omega$, and we have to move it into the
non-linear part.
This is taken care of in the evolver

```python
bec.evolve_comoving_dGPE(self, number_of_steps, velx,method='ETD2RK')
```

Here it is assumed that the boost is in the $x$-direction, and that the dissipative factor is spatially dependent.

[^andersonObservationBoseEinsteinCondensation1995]: Anderson, M. H., Ensher, J. R., Matthews, M. R., Wieman, C. E., & Cornell, E. A. (1995). Observation of Bose-Einstein Condensation in a Dilute Atomic Vapor. Science, 269(5221), 198–201. [https://doi.org/10.1126/science.269.5221.198](https://doi.org/10.1126/science.269.5221.198)
[^dalfovo1999theory]: Dalfovo, F., Giorgini, S., Pitaevskii, L. P. and Stringari, S. (1999). Theory of Bose-Einstein condensation in trapped gases. Reviews of Modern Physics. 71, 3, 463. [https://doi.org/10.1103/RevModPhys.71.463](https://doi.org/10.1103/RevModPhys.71.463)
[^kevrekidis2007emergent]: Kevrekidis, P. G.,  Frantzeskakis, D. J. and  Carretero-González, R. (2008). Emergent nonlinear phenomena in Bose-Einstein condensates: theory and experiment. Springer Science & Business Media. Berlin.
[^pitaevskiiBook]: Pitaevskii, L. and Stringari, S. (2016). Bose-Einstein Condensation and Superfluidity. Oxford University Press. [https://doi.org/10.1093/acprof:oso/9780198758884.001.0001](https://doi.org/10.1093/acprof:oso/9780198758884.001.0001})
[^gardiner2003stochastic]: Gardiner, C. W. and Davis, M. J. (2003). The stochastic Gross-Pitaevskii equation: II. Journal of Physics B: Atomic, Molecular and Optical Physics. 36, 23, 4731. [https://doi.org/10.1088/0953-4075/36/23/010](https://doi.org/10.1088/0953-4075/36/23/010)
[^rooney2012stochastic]: Rooney, S. J., Blakie, P. B. and Bradley, A. S. (2012). Stochastic projected Gross-Pitaevskii equation. Physical Review A. 86, 5, 053634. [https://doi.org/10.1103/PhysRevA.86.053634](https://doi.org/10.1103/PhysRevA.86.053634)
[^bradley2012energy]: Bradley, A. S. and Anderson, B. P. (2012). Energy spectra of vortex distributions in two-dimensional quantum turbulence. Physical Review X. 2, 4, 041001 [https://doi.org/10.1103/PhysRevX.2.041001](https://doi.org/10.1103/PhysRevX.2.041001)
[^skaugenUnifiedPerspectiveTwodimensional2018]: Skaugen, A. (2018). A Unified Perspective on Two-Dimensional Quantum Turbulence and Plasticity. PhD Thesis, University of Oslo. [http://urn.nb.no/URN:NBN:no-69394](http://urn.nb.no/URN:NBN:no-69394)
[^minguzzi2004numerical]: Minguzzi, A., Succi, S., Toschi, F., Tosi, M. P. and Vignolo, P. (2004). Numerical methods for atomic quantum gases with applications to Bose-Einstein condensates and to ultracold fermions. Physics reports. 395, 4-5, 223-355. [https://doi.org/10.1016/j.physrep.2004.02.001](https://doi.org/10.1016/j.physrep.2004.02.001)
[^nore1997kolmogorov]: Nore, C., Abid, M. and Brachet, M. E. (1997). Kolmogorov turbulence in low-temperature superflows. Physical review letters. 78, 20, 3896. [https://doi.org/10.1103/PhysRevLett.78.3896](https://doi.org/10.1103/PhysRevLett.78.3896)
[^ronning2020classical]: Rønning, J., Skaugen, A., Hernández-García, E., López, C. and Angheluta, L. (2020). Classical analogies for the force acting on an impurity in a Bose-Einstein condensate. New Journal of Physics. 22, 7, 073018. [https://doi.org/10.1088/1367-2630/ab95de](https://doi.org/10.1088/1367-2630/ab95de)
[^astrakharchik2004motion]: Astrakharchik, G. E. and Pitaevskii, L. P. (2004). Motion of a heavy impurity through a {Bose-Einstein} condensate. Physical Review A. 70, 1, 013608. [https://doi.org/10.1103/PhysRevA.70.013608](https://doi.org/10.1103/PhysRevA.70.013608)
[^pinsker2017gaussian]: Pinsker, F. (2017). Gaussian impurity moving through a {Bose-Einstein} superfluid. Physica B: Condensed Matter. 521, 36-42. [https://doi.org/10.1016/j.physb.2017.06.038](https://doi.org/10.1016/j.physb.2017.06.038)
[^Pismen]: Pismen, L.M. (1999). Vortices in Nonlinear Fields: From Liquid Crystals to Superfluids, From Non-Equilibrium Patterns to Cosmic Strings. Oxford university press. Oxford
[^reeves2015identifying]: Reeves, M. T., Billam, T. P., Anderson, B. P. and Bradley, A. S. (2015). Identifying a superfluid Reynolds number via dynamical similarity. Physical Review Letters. 114, 15, 155302. [https://doi.org/10.1103/PhysRevLett.114.155302](https://doi.org/10.1103/PhysRevLett.114.155302)
