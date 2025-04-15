# Class: Quantum Mechanics

Quantum mechanics describes the behavior of nature at the scale of atoms and subatomic particles. The Schrödinger equation is a fundamental partial differential equation that governs how the quantum state of a physical system evolves over time.

This class simulates quantum mechanics by evolving the Schrödinger equation.

```python
file: comfit/quantum_mechanics/quantum_mechanics.py 
class: QuantumMechanics
```

See the ComFiT Library Reference below for a complete list of class methods and their usage.

<div class="grid cards" style="display: flex; flex-wrap: wrap;">
    <a href="https://comfitlib.com/library_reference/quantum_mechanics/" class="card" style="min-width: 160px; flex: 0 1 calc(100.00% - 10px); margin: 5px;">
        <div style="text-align: center;">
            <strong> ComFiT Library Reference </strong>
        </div>
    </a>
</div>

## Example

The following example demonstrates how to set up a 1D quantum system with a Gaussian wave packet and a potential barrier. It runs smoothly with `comfit 1.8.4`.

```python
import comfit as cf
import matplotlib.pyplot as plt
import numpy as np

# Set up a 1D quantum system
qm = cf.QuantumMechanics(1, xlim=[-50,50], xRes=1001, dt=0.1)

# Initialize a Gaussian wavepacket at x=5 with velocity=1
qm.conf_initial_condition_Gaussian(position=5, width=1, initial_velocity=1)

# Add a potential barrier (optional)
qm.V_ext = 0.5 * (qm.x > 10) * (qm.x < 12)  # Barrier from x=10 to 12

height = np.max(abs(qm.psi))  # Get the maximum value of the wavefunction for plotting limits

# Optional: Animate it
for n in range(61):
    qm.evolve_schrodinger(5)
    fig, ax = qm.plot_complex_field(qm.psi)
    qm.plot_field(qm.V_ext, fig=fig, ax=ax, ylim=[0,height], xlim=[0,20])
    qm.plot_save(fig, n)
cf.tool_make_animation_gif(n)  # Creates animation.gif
```

![Quantum Mechanics](../img/quantum_mechanics_barrier_reflection.gif#only-light)
![Quantum Mechanics](../img/quantum_mechanics_barrier_reflection-colorinverted.gif#only-dark)

## The Schrödinger equation

The time-dependent Schrödinger equation for a particle of mass $m_e$ in a potential $V$ is given by:

$$
\mathfrak i \hbar \partial_t \psi = \left [-\frac{\hbar^2}{2m_e} \nabla^2 + V \right ] \psi.
$$

where $\psi$ is the wave function, $\hbar$ is the reduced Planck constant.

To simplify the equation for numerical simulation, we can introduce dimensionless units. Dividing by the Hartree energy $E_h = \frac{\hbar^2}{m_e a_0^2}$ (where $a_0$ is the Bohr radius) yields:

$$
\mathfrak i \frac{\hbar}{E_h} \partial_t \psi = \frac{1}{E_h}\left [-\frac{\hbar^2}{2m_e} \nabla^2 + V \right ] \psi.
$$

Expressing time in units of $\tau = \frac{\hbar}{E_h}$, potential energy $V$ in units of $E_h$, and length squared in units of $a_0^2$, we obtain the dimensionless Schrödinger equation:

$$
\partial_t \psi = \mathfrak i\left [\frac{1}{2} \nabla^2 - V \right ] \psi.
$$

This is the form implemented in the `QuantumMechanics` class. The linear operator $\omega$ and its Fourier space representation $\omega_{\mathfrak f}$ are:

$$
\omega = \mathfrak i\frac{1}{2} \nabla^2 \quad \Rightarrow \quad \omega_{\mathfrak f}(\mathbf{k}) = -\mathfrak i \frac{1}{2} \mathbf k^2
$$

The non-linear part (in this context, the potential term which depends on $\psi$ only through multiplication) is $N = -\mathfrak i V \psi$.

| Atomic unit of | Value               |
|----------------|---------------------|
| Length ($a_0$) | 0.529 Å (Angstrom)  |
| Energy ($E_h$) | 27.2 eV (electron volts)  |
| Time ($\tau$)  | 24.2 as (attoseconds)    |

## The Born rule

The Born rule states that the probability density of finding the particle at position $\mathbf{r}$ at time $t$ is given by $|\psi(\mathbf{r}, t)|^2$. Consequently, the probability $P$ of finding the particle within a volume $\mathcal{V}$ is:

$$
P = \int_{\mathcal{V}} d^d r |\psi(\mathbf{r}, t)|^2.
$$

The total probability of finding the particle anywhere in space must be 1, leading to the normalization condition:

$$
\int d^d r |\psi(\mathbf{r}, t)|^2 = 1.
$$

## The Momentum representation

In quantum mechanics, the Fourier transform provides the connection between the position representation ($\psi(\mathbf{r})$) and the momentum representation ($\phi(\mathbf{k})$) of a quantum state. This relationship highlights the wave-particle duality inherent in quantum theory.

The momentum-space wave function $\phi(\mathbf{k})$, where $|\phi(\mathbf{k})|^2$ represents the probability density of finding the particle with wave vector $\mathbf{k}$ (momentum $\mathbf{p} = \hbar \mathbf{k}$), is related to the position-space wave function $\psi(\mathbf{r})$ via the Fourier transform:

$$\phi(\mathbf{k}) = \sqrt{(2\pi)^d} \, \psi_{\mathfrak f} (\mathbf{k})$$

where $\psi_{\mathfrak f} (\mathbf{k})$ is the Fourier transform of $\psi(\mathbf{r})$ as defined in the [Base System documentation](https://comfitlib.com/ClassBaseSystem/#fourier-transformations), and $d$ is the spatial dimension. The factor $\sqrt{(2\pi)^d}$ ensures that if $\psi(\mathbf{r})$ is normalized in position space, then $\phi(\mathbf{k})$ is normalized in momentum space:

$$
\int d^d k |\phi(\mathbf{k})|^2 = 1.
$$

### Physical Significance of the Fourier Transform in QM

1.  **Complementarity**: The Fourier transform relationship mathematically embodies Heisenberg's uncertainty principle. A wave function highly localized in position space (narrow $\psi(\mathbf{r})$) corresponds to a widely spread momentum distribution (broad $\phi(\mathbf{k})$), and vice versa.
2.  **Operator Correspondence**: In the position representation, the momentum operator is $\hat{\mathbf{p}} = -\mathfrak i\hbar\nabla$. The Fourier transform maps this differential operator in position space to a simple multiplicative operator ($\hbar \mathbf{k}$) in momentum space.
3.  **Energy Eigenstates**: For a free particle ($V=0$), the energy eigenstates are plane waves, $\psi(\mathbf r) \propto e^{\mathfrak i \mathbf k \cdot \mathbf r}$, which are also momentum eigenstates.

## A wave packet

A Gaussian wave function, often referred to as a wave packet, provides a localized representation of a particle in quantum mechanics, balancing position and momentum uncertainty. A normalized Gaussian wave packet centered at $\mathbf{r}_0$ with width $\sigma$ is given by:

$$
\psi(\mathbf r) = \sqrt{( 2\pi \sigma )^{-d/2} \exp\left ({-\frac{(\mathbf r - \mathbf r_0)^2} {(2\sigma^2)}}\right )} ,
$$

such that the probability density $|\psi(\mathbf{r})|^2 = (2\pi \sigma^2)^{-d/2} \exp\left ({-\frac{(\mathbf r - \mathbf r_0)^2} {2\sigma^2}}\right )$ is a Gaussian distribution with standard deviation $\sigma$.

An initial average velocity $\mathbf v_0$ (corresponding to an average momentum $\hbar \mathbf{k}_0 = m_e \mathbf{v}_0$) can be imparted to the wave packet by multiplying it with a complex phase factor $e^{\mathfrak i \mathbf{k}_0 \cdot \mathbf r} = e^{\mathfrak i (m_e/\hbar)\mathbf v_0 \cdot \mathbf r}$. In dimensionless units where $m_e=\hbar=1$, this simplifies to $e^{\mathfrak i \mathbf v_0 \cdot \mathbf r}$.

Such an initial condition can be configured using the function `qm.conf_initial_condition_Gaussian`.
