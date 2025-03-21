# Class: Quantum Mechanics

Quantum mechanics is one of the most classic examples of field theories in physics.
The Schrödinger equation is a partial differential equation that describes how the quantum state of a physical system changes over time.

In this class, we simulate quantum mechanics through the evolution of the Schrödinger equation.

```python
file: comfit/models/quantum_mechanics.py 
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

The following example demonstrates how to set up a 1D quantum system with a Gaussian wave packet and a potential barrier.
It runs smoothly with `comfit 1.8.4`.

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

height = np.max(abs(qm.psi))  # Get the maximum value of the wavefunction

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

Evolving according to the Schrödinger equation with electron mass

$$
\mathfrak i \hbar \partial_t \psi = \left [-\frac{\hbar^2}{2m_e} \nabla^2 + V \right ] \psi.
$$

Dividing by the Hartree energy $E_h = \frac{\hbar^2}{m_e a_0^2}$

$$
\mathfrak i \frac{\hbar}{E_h} \partial_t \psi = \frac{1}{E_h}\left [-\frac{\hbar^2}{2m_e} \nabla^2 + V \right ] \psi.
$$

When expressing time in units of $\frac{\hbar}{E_h}$, potential energy in units of $E_h$ and length squared in units of 

$$
\frac{\hbar^2}{E_h m_e} = \frac{\hbar^2}{\frac{\hbar^2}{m_e a_0^2} m_e} = a_0^2,
$$

we get the Schrödinger equation in its dimensionless form

$$
 \partial_t \psi = \mathfrak i\left [\frac{1}{2} \nabla^2 - V \right ] \psi.
$$

So 

$$
\omega = \mathfrak i\frac{1}{2} \nabla^2
\Rightarrow \omega_{\mathfrak f} = -\mathfrak i \frac{1}{2} \mathbf k^2
$$

| Atomic unit of | Value               |
|----------------|---------------------|
| Length         | 0.529 Å (Angstrom)  |
| Energy         | 27.2 eV (electron volts)  |
| Time           | 24.2 aS (atto seconds)    |


## The Born rule

The Born rule states that the probability $p$ of measuring a particle in the interval $[a,b]$ is given by 

$$
p = \int_a^b dx |\psi(x)|^2.
$$

## The Momentum representation

In quantum mechanics, the Fourier transform serves as the mathematical bridge between the position and momentum representations of a quantum state. 
This relationship reveals the wave-particle duality at the heart of quantum theory.

The same quantum state can alternatively be described in momentum space by the wave function $\psi(k)$, where $|\phi(k)|^2$, where $k$ is the wave number, gives the probability density of finding the particle with momentum $\hbar k$.
The momentum wave function $\phi(p)$ is related to the position wave function $\psi(x)$ through the Fourier transform:

$$\phi(k) = \sqrt{2\pi} \psi_{\mathfrak f} (k)$$

with the Fourier transform $\psi_{\mathfrak f} (k)$ defined as in [The Base System documentation](https://comfitlib.com/ClassBaseSystem/#fourier-transformations).
The extra factor of $\sqrt{2\pi}$ is a comes from our definition of the Fourier transform, see below.

### Normalization of the Fourier Transform

Consider the Fourier transform in $n$-dimensions as defined:

$$\psi_{\mathfrak f}(\mathbf{k}) = \mathcal F[\psi] = \frac{1}{(2\pi)^n} \int d^n r\, e^{-i\mathbf{k}\cdot\mathbf{r}} \psi(\mathbf{r}),$$

with the inverse:

$$\psi(\mathbf{r}) = \mathcal F^{-1}[\psi_{\mathfrak f}] = \int d^n k\, e^{i\mathbf{k}\cdot\mathbf{r}} \psi_{\mathfrak f}(\mathbf{k}),$$

where $d^n r$ and $d^n k$ denote integration over $n$-dimensional position and wave vector spaces, respectively. In quantum mechanics, the position-space wave function $\psi(\mathbf{r})$ is normalized such that:

$$\int d^n r\, |\psi(\mathbf{r})|^2 = 1,$$

ensuring the total probability is 1 across the $n$-dimensional space.

For the Fourier representation, substituting the inverse transform into the normalization condition and evaluating the inner integral yields a Dirac delta, $\int d^n r\, e^{i(\mathbf{k}-\mathbf{k'})\cdot\mathbf{r}} = (2\pi)^n \delta^n(\mathbf{k}-\mathbf{k'})$. This leads to:

$$\int d^n k\, (2\pi)^n |\psi_{\mathfrak f}(\mathbf{k})|^2 = 1,$$

or:

$$\int d^n k\, |\psi_{\mathfrak f}(\mathbf{k})|^2 = \frac{1}{(2\pi)^n}.$$

The factor $(2\pi)^n$ arises from the convention placing $\frac{1}{(2\pi)^n}$ in the forward transform.

To define a $\mathbf{k}$-space wave function $\phi(\mathbf{k})$ where $|\phi(\mathbf{k})|^2$ is the probability density per unit $\mathbf{k}$ in $n$-dimensions, satisfying:

$$\int d^n k\, |\phi(\mathbf{k})|^2 = 1,$$

set:

$$\phi(\mathbf{k}) = \sqrt{(2\pi)^n} \, \psi_{\mathfrak f}(\mathbf{k}).$$

Then:

$$|\phi(\mathbf{k})|^2 = (2\pi)^n |\psi_{\mathfrak f}(\mathbf{k})|^2,$$

and:

$$\int d^n k\, |\phi(\mathbf{k})|^2 = (2\pi)^n \int d^n k\, |\psi_{\mathfrak f}(\mathbf{k})|^2 = (2\pi)^n \cdot \frac{1}{(2\pi)^n} = 1.$$

Thus, $|\phi(\mathbf{k})|^2$ serves as the probability density in $\mathbf{k}$-space, adjusted for this Fourier transform convention.
In momentum space ($\mathbf{p} = \hbar \mathbf{k}$), additional scaling by $\hbar^n$ may apply, but $\phi(\mathbf{k}) = \sqrt{(2\pi)^n}  \psi_{\mathfrak f}(\mathbf{k})$ ensures proper normalization in $\mathbf{k}$-space.

### Physical Significance

1. **Complementarity**: The Fourier transform relationship embodies Heisenberg's uncertainty principle.
The more localized a wave function is in position space, the more spread out its Fourier transform is in momentum space, and vice versa.

2. **Operator Correspondence**: In the position representation, the momentum operator is $\hat{p} = -\mathfrak i\hbar\frac{d}{dx}$.
The Fourier transform converts this differential operator to a multiplication operator in momentum space.

3. **Energy Eigenstates**: For a free particle, the energy eigenstates are momentum eigenstates, which are plane waves in position space: $\psi(x) \propto e^{\mathfrak ikx/\hbar}$.


## A wave packet

A Gaussian wave function is often called a wave packet and can visualize the position and motion of a particle in a quantum mechanical system.
An initial wave function with a widt of $\sigma$ is given by

$$
\psi(\mathbf r) = \sqrt{( 2\pi \sigma )^{-d/2} \exp\left ({-\frac{(\mathbf r - \mathbf r_0)^2} {(2\sigma^2)}}\right )} ,
$$

so that $|\psi|^2$ is a Gaussian distribution.
An initial velocity $\mathbf v_0$ can be given to the wave packet by multiplying with a complex phase $e^{\mathfrak i \mathbf v_0 \cdot \mathbf r}$.
Such an initial condition can be configured by the function `qm.conf_initial_condition_Gaussian`.