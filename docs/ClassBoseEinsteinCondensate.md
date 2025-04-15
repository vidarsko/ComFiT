# Class: Bose-Einstein Condensate

There are two types of particles in the world: fermions and bosons. Whereas fermions can never occupy the same quantum state due to the Pauli Exclusion Principle, the same is not true for bosons. A Bose-Einstein condensate (BEC) is a state of matter consisting of ultra-cold bosons that undergo a phase transition at a low critical temperature, causing most bosons to occupy the lowest quantum energy state (the ground state) of the system. It was theorized by Satyendra Nath Bose and Albert Einstein in the 1920s as a new state of matter and was first produced experimentally in 1995 by Eric Cornell and Carl Wieman[^andersonObservationBoseEinsteinCondensation1995].

This class simulates a Bose-Einstein condensate in 1, 2, and 3 dimensions using the Gross-Pitaevskii equation (GPE).

```python
file: comfit/models/bose_einstein_condensate.py 
class: BoseEinsteinCondensate
```

See the ComFiT Library Reference below for a complete list of class methods and their usage.

<div class="grid cards" style="display: flex; flex-wrap: wrap;">
    <a href="https://comfitlib.com/library_reference/bose_einstein_condensate/" class="card" style="min-width: 160px; flex: 0 1 calc(100.00% - 10px); margin: 5px;">
        <div style="text-align: center;">
            <strong> ComFiT Library Reference </strong>
        </div>
    </a>
</div>

## Variables and Parameters

The primary field in the Bose-Einstein condensate model is the complex wave function $\psi$:

```python
bec.psi
```

The `BoseEinsteinCondensate` class accepts the same keyword arguments as the `BaseSystem` class, with the addition of the following specific parameter:

| Keyword | Definition         | Default Value |
|---------|--------------------|---------------|
| `gamma` | Dissipative factor | $0.01$        |

This parameter allows for the customization and fine-tuning of the Bose-Einstein condensate simulation.

## Model

The BEC in the mean-field regime is described by the Gross-Pitaevskii equation (GPE)[^dalfovo1999theory] [^kevrekidis2007emergent]. This is a non-linear Schrödinger equation which reads:

$$
i\hbar \partial_t\psi = \left[-\frac{\hbar^2}{2m} \nabla^2+ V_{ext}(\mathbf{r}, t) -\mu +g|\psi|^2 \right]\psi.
$$

Here, $\mu$ is the chemical potential, $m$ is the mass of the bosons, $g$ is an interaction parameter, and $\hbar$ is the reduced Planck constant. $\psi$ is the wave function describing the condensate phase, and $V_{ext}(\mathbf{r}, t)$ is an external potential, which can depend on position $\mathbf{r}$ and time $t$. The GPE can be obtained from the variational principle $\mathfrak i \hbar \partial_t \psi = \frac{\delta K}{\delta \psi^*}$ [^kevrekidis2007emergent][^pitaevskiiBook] with the Hamiltonian:

$$
K = \int d\mathbf{r} \left[\frac{\hbar^2}{2m}|\nabla\psi|^2 +(V_{ext}(\mathbf{r}, t) -\mu)|\psi|^2 +\frac{g}{2}|\psi|^4 \right].
$$

We introduce dimensionless units for length $\xi = \hbar/\sqrt{m\mu}$ and time $\tau = \xi/c$, and rescale the wave function as $\psi \rightarrow \sqrt{\frac{g}{\mu}}\psi$. Additionally, we include a dissipative factor $\gamma$. This results in the damped Gross-Pitaevskii equation (dGPE) in dimensionless form as [^gardiner2003stochastic][^rooney2012stochastic][^bradley2012energy][^skaugenUnifiedPerspectiveTwodimensional2018]:

!!! equation "The damped Gross-Pitaevski equation (`evolve_dGPE`)"
    $$
    \partial_t \psi = -(\mathfrak{i}+\gamma) \left[-\frac{1}{2}\nabla^2 + V_{ext}(\mathbf{r}, t) -1 +|\psi|^2 \right]\psi.
    $$

We can split this equation into its linear and non-linear parts:

$$
\omega = (\mathfrak{i}+\gamma) (1 + \frac{1}{2}\nabla^2), \quad N = - (\mathfrak{i} + \gamma) (V_{ext}(\mathbf{r}, t) + |\psi|^2)\psi
$$

The Hamiltonian and its density in dimensionless units are calculated by the functions:

```python
bec.calc_hamiltonian()
bec.calc_hamiltonian_density()
```

## Approximation of Ground States

In simulations, it is often convenient to start in a configuration that is close to the ground state. We can estimate this ground state by noting that the GPE dissipates energy when evolved in imaginary time ($t \rightarrow -it$) [^minguzzi2004numerical][^kevrekidis2007emergent]. We refer to this as *evolving the dGPE in imaginary time*. Given an external potential $V_{ext}$, we can therefore find an approximation to the ground state by starting with an initial guess and then removing energy from the guessed state by evolving the equations in imaginary time.

To obtain an initial guess for the ground state for a given potential $V_{ext}$, we use the Thomas-Fermi approximation [^dalfovo1999theory][^kevrekidis2007emergent]. In this approximation, we assume that $\psi$ is slowly varying, allowing us to neglect the Laplacian term ($\nabla^2 \psi \approx 0$). Looking for stationary solutions ($\partial_t \psi = 0$) to the dGPE, we obtain the equation:

$$
0 = (V_{ext}(\mathbf{r}) - 1 +|\psi|^2 )\psi.
$$

This equation has two solutions: $\psi = 0$ and $|\psi|^2 = 1-V_{ext}(\mathbf{r})$. In the case where $V_{ext}(\mathbf{r}) > 1$, only $\psi = 0$ is physically possible (since density $|\psi|^2$ must be non-negative). In the case where $V_{ext}(\mathbf{r}) < 1$, both solutions exist. To determine which represents the ground state, we consider the Hamiltonian in the Thomas-Fermi limit:

$$
K_{TF} = \int d\mathbf{r} \left[(V_{ext}(\mathbf{r}) -1)|\psi|^2 + \frac{1}{2} |\psi|^4 \right].
$$

Inserting the ansatz $\psi = 0$ yields $K_{TF} = 0$. Inserting the ansatz $|\psi|^2 = 1-V_{ext}(\mathbf{r})$ gives $K_{TF} = \int d\mathbf{r} [-(1-V_{ext})^2 + \frac{1}{2}(1-V_{ext})^2] = -\frac{1}{2} \int d\mathbf{r} (1-V_{ext})^2$, which is negative for $V_{ext}(\mathbf{r}) < 1$. We therefore conclude that the Thomas-Fermi ground state approximation is:

$$
\psi_{TF}(\mathbf{r}) = \begin{cases}
    0 & \text{if } V_{ext}(\mathbf{r}) \ge 1 \\
    \sqrt{1 -V_{ext}(\mathbf{r})}  & \text{if } V_{ext}(\mathbf{r}) < 1
    \end{cases}
$$

This ground state guess can be initialized using:

```python
bec.conf_initial_condition_Thomas_Fermi()
```

Once this guess is initialized, the wave function can be propagated in imaginary time to relax towards the true ground state using the function:

```python
bec.evolve_relax(number_of_steps, method='ETD2RK')
```

Note that the external potential $V_{ext}$ must be constant during this relaxation evolution.

## Potentials

Many interesting dynamics in a BEC arise from its interaction with an external potential. This potential is represented by the function:

```python
bec.V_ext(t) 
```

By default, this is initialized as a time-independent potential returning 0 everywhere. The potential can be changed using the function:

```python
bec.conf_external_potential(V_ext, additive=False)
```

This function can set the potential based on whether `V_ext` is a callable function (expecting time `t` as an argument), a constant float, or a NumPy array representing a static potential field. If `additive=True`, the provided `V_ext` (if constant or an array) is added to the existing potential. If `V_ext` is a function, it should have the signature:

```python
def V(t):
     # calculate potential based on time t
     # requires access to spatial coordinates like bec.x, bec.y, bec.z
     potential_field = ... 
     return potential_field
```

The evolver will evaluate this function using the current simulation time `bec.time`. An example using a Gaussian stirring potential is provided in the example folder.

To simplify setup, some common potentials are provided:

* **Harmonic Potential:** Configured using `bec.conf_harmonic_potential(trapping_strength)`. The potential takes the form $V_H(\mathbf{r}) = \texttt{trapping_strength} \times |\mathbf{r} - \mathbf{r}_{mid}|^2$.
* **Gaussian Potential:** Can be created using `bec.calc_Gaussian(position, width, top=strength)`. This returns $V_g(\mathbf{r}) = \texttt{strength} \times e^{-|\mathbf{r} - \mathbf{r}_{position}|^2/(2 \times \texttt{width}^2)}$. You would then use `conf_external_potential` to set this as the external potential.

Much of the interesting physics occurs when the potential is time-dependent. For this, one can define a function `V(t)` as shown above and update the external potential by calling `bec.conf_external_potential(V)`.

## Hydrodynamics

The dGPE can be transformed into a hydrodynamic description of the BEC. The first step is to introduce the Madelung transformation $\psi = \sqrt{\rho} e^{i\theta}$, where $\rho = |\psi|^2$ is the superfluid density [^kevrekidis2007emergent]. For $\gamma = 0$, this density is conserved and satisfies the continuity equation:

$$
\partial_t \rho + \nabla\cdot \mathbf{J}_s = 0,
$$

with the superfluid current given by:

$$
\mathbf{J}_s = \Im(\psi^* \nabla \psi) = \rho \nabla \theta = \rho \mathbf{v}_s.
$$

Here, the superfluid velocity is introduced as $\mathbf{v}_s = \nabla \theta$. The superfluid current can be calculated using the function:

```python
bec.calc_superfluid_current()   
```

We can also substitute the Madelung transformation into the Hamiltonian to get [^bradley2012energy][^nore1997kolmogorov]:

$$
K = \int d\mathbf{r} \left[\frac{1}{2}\rho v_s^2 +\frac{1}{8} \frac{|\nabla \rho|^2}{\rho} + (V_{ext}-1)\rho +\frac{1}{2}\rho^2 \right].
$$
Note: The original text had $\rho^4$, but based on the GPE Hamiltonian, it should be $\rho^2$ (corresponding to $|\psi|^4$).

The first term, $\frac{1}{2}\rho v_s^2$, represents the kinetic energy density of the condensate. To calculate the total kinetic energy, it is convenient to introduce the density-weighted velocity $\mathbf{u} = \sqrt{\rho}\mathbf{v}_s$ [^bradley2012energy]. This quantity has the advantage of not being singular at the core of topological defects (where $\rho \to 0$). Using this, the kinetic energy can be written as $E_k = \int d\mathbf{r} \frac{1}{2}|\mathbf{u}|^2$. The density-weighted velocity and the total kinetic energy can be calculated using the functions:

```python
bec.calc_velocity()
bec.calc_kinetic_energy()
```

Furthermore, by inserting the Madelung transformation into the dGPE, one can derive equations resembling the Navier-Stokes equations [^kevrekidis2007emergent][^bradley2012energy]:

$$\begin{aligned}
    \partial_t \rho + \nabla\cdot(\rho \mathbf{v}_s) &= 2\gamma \rho (1-V_{eff}), \\
    \partial_t \mathbf{v}_s + (\mathbf{v}_s \cdot \nabla) \mathbf{v}_s &=  - \nabla(V_{ext} + \rho - \frac{1}{2\sqrt{\rho}}\nabla^2\sqrt{\rho}) + \frac{\gamma}{2}\nabla^2 \mathbf{v}_s.
\end{aligned}$$

Note: The exact form of the momentum equation can vary depending on the derivation details and approximations (e.g., quantum pressure term $\frac{1}{2\sqrt{\rho}}\nabla^2\sqrt{\rho}$). The density equation shows that the condensate density is only conserved when $\gamma = 0$. $V_{eff}$ includes the external potential and interaction terms.

## Forces on External Potential

To study how the BEC interacts with impurities, one can model the impurity as an external potential (e.g., Gaussian) and measure the force exerted on it [^ronning2020classical][^astrakharchik2004motion][^pinsker2017gaussian]. According to the Ehrenfest theorem, the force exerted by the external potential on the condensate is given by $\mathbf{F}_{cond} = -\langle \nabla V_{ext}\rangle$. By Newton's third law, the force exerted by the condensate on the potential (or the object creating it) is $\mathbf{F}_{pot} = - \mathbf{F}_{cond} = \langle \nabla V_{ext}\rangle$. Written explicitly, this is:

$$
\mathbf{F}_{pot} = \int d\mathbf{r} |\psi|^2 \nabla V_{ext}(\mathbf{r}) = -\int d\mathbf{r} V_{ext}(\mathbf{r}) \nabla(|\psi|^2),
$$

where the second equality follows from integration by parts (assuming boundary terms vanish). This force is calculated by the function:

```python
bec.calc_force_on_external_potential()
```

Note that this calculates the total force exerted by the condensate on the source of the *entire* external potential $V_{ext}$.

## Quantized Vortices

Topological defects in a BEC manifest as quantized vortices. This arises because the superfluid velocity is the gradient of the phase, $\mathbf{v}_s = \nabla \theta$. Consequently, the circulation around any closed loop $C$ is quantized:

$$
\Gamma = \oint_C d\mathbf{l} \cdot \mathbf{v}_s = \oint_C d\mathbf{l} \cdot \nabla \theta = \Delta\theta = 2\pi n,
$$

where $n$ must be an integer for the wave function $\psi = \sqrt{\rho}e^{i\theta}$ to be single-valued. In two dimensions, these vortices are point defects, while in three dimensions, they are line defects. The locations (nodes) of these vortices can be identified using the function:

```python
bec.calc_vortex_nodes(dt_psi=None)
```

This function finds the vortex nodes in two or three dimensions. If the time derivative of the wave function, $\partial_t \psi$, is provided as the argument `dt_psi`, the function also calculates the velocity of the defects.

Vortices can be created in several ways: by manually inserting them into the initial state (using functions like `conf_insert_vortex`, `conf_insert_vortex_dipole`, `conf_insert_vortex_ring`), by stirring the condensate with a moving potential, or by relaxing a disordered initial state (quenching).

## Comoving Frame

When studying a BEC stirred by a potential moving at velocity $-\mathbf{V}_p$, it is sometimes convenient to switch to a reference frame moving with the potential (the comoving frame). In this frame, the potential is stationary, and the BEC flows past it with velocity $\mathbf{V}_p$. The transformed dGPE reads:

$$
\partial_t \psi = -\mathbf{V}_p \cdot \nabla \psi - (\mathfrak{i}+\gamma) \left[-\frac{1}{2}\nabla^2 + V_{ext}(\mathbf{r}) -1 +|\psi|^2 \right]\psi.
$$

Note that a standard Galilean transformation of the GPE ($\gamma=0$) often includes a phase factor $e^{i(\mathbf{V}_p \cdot \mathbf{r} - \frac{1}{2} V_p^2 t)}$ to ensure the transformed velocity represents the velocity in the new frame[^Pismen]. However, the dGPE ($\gamma \neq 0$) is not Galilean invariant, so we omit this phase factor here. Consequently, the superfluid velocity calculated from $\psi$ in this frame still represents the velocity in the original lab frame [^ronning2020classical].

To prevent artificial recycling of excitations from the outflow boundary back into the incoming flow due to periodic boundary conditions, a dissipative buffer region can be introduced around the computational domain where $\gamma$ is significantly larger than in the bulk [^reeves2015identifying]. The dissipative factor becomes a function of space, typically defined as:

$$
\gamma(\mathbf{r}) = \max[\gamma_x(x), \gamma_y(y), \gamma_z(z)] \quad \text{(3D)}
$$

or

$$
\gamma(\mathbf{r}) = \max[\gamma_x(x), \gamma_y(y)] \quad \text{(2D)}
$$

where, for example:

$$
\gamma_x(x) = \gamma_0 + \frac{\gamma_{max}-\gamma_0}{2} \left( \tanh\left(\frac{x - (x_{mid} + w_x/2)}{d}\right) - \tanh\left(\frac{x - (x_{mid} - w_x/2)}{d}\right) + 2 \right).
$$

Here, $\gamma_0$ is the bulk dissipation value, $\gamma_{max}$ is the high value in the buffer (often implicitly set, e.g., $\gamma_{max} \approx \gamma_0 + 1$), $d$ is the width of the interface between the bulk and the buffer, and $w_x$ defines the size of the central low-dissipation region along the x-axis. Similar functions define $\gamma_y(y)$ and $\gamma_z(z)$. When $\gamma$ depends on space, the term involving $\gamma \nabla^2 \psi$ can no longer be treated as part of the linear operator $\omega$ and must be handled within the non-linear term calculation. This is managed by the specific evolver:

```python
bec.evolve_comoving_dGPE(number_of_steps, velx, method='ETD2RK')
```

This evolver assumes the boost (and thus the flow) is in the positive x-direction with speed `velx` and correctly handles the spatially dependent dissipation $\gamma(\mathbf{r})$.

[^andersonObservationBoseEinsteinCondensation1995]: Anderson, M. H., Ensher, J. R., Matthews, M. R., Wieman, C. E., & Cornell, E. A. (1995). Observation of Bose-Einstein Condensation in a Dilute Atomic Vapor. Science, 269(5221), 198–201. [https://doi.org/10.1126/science.269.5221.198](https://doi.org/10.1126/science.269.5221.198)
[^dalfovo1999theory]: Dalfovo, F., Giorgini, S., Pitaevskii, L. P. and Stringari, S. (1999). Theory of Bose-Einstein condensation in trapped gases. Reviews of Modern Physics. 71, 3, 463. [https://doi.org/10.1103/RevModPhys.71.463](https://doi.org/10.1103/RevModPhys.71.463)
[^kevrekidis2007emergent]: Kevrekidis, P. G., Frantzeskakis, D. J. and Carretero-González, R. (2008). Emergent nonlinear phenomena in Bose-Einstein condensates: theory and experiment. Springer Science & Business Media. Berlin.
[^pitaevskiiBook]: Pitaevskii, L. and Stringari, S. (2016). Bose-Einstein Condensation and Superfluidity. Oxford University Press. [https://doi.org/10.1093/acprof:oso/9780198758884.001.0001](https://doi.org/10.1093/acprof:oso/9780198758884.001.0001)
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
