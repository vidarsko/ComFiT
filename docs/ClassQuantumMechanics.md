# Class: Quantum Mechanics

Quantum mechanics is one of the most classic examples of field theories in physics.
The Schrödinger equation is a partial differential equation that describes how the quantum state of a physical system changes over time.

In this class, we simulate quantum mechanics through the evolution of the Schrödinger equation.

```python
file: comfit/models/quantum_mechanics.py 
class: QuantumMechanics
```

## Functions

Let `qm` be an instance of the class `QuantumMechanics`.

- `qm.evolve_schrodinger()`: 

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

## A wave packet

A Gaussian wave function is often called a wave packet and can visualize the position and motion of a particle in a quantum mechanical system.
An initial wave function with a widt of $\sigma$ is given by

$$
\psi(\mathbf r) = \sqrt{( 2\pi \sigma )^{-d/2} \exp\left ({-\frac{(\mathbf r - \mathbf r_0)^2} {(2\sigma^2)}}\right )} ,
$$

so that $|\psi|^2$ is a Gaussian distribution.
An initial velocity $\mathbf v_0$ can be given to the wave packet by multiplying with a complex phase $e^{\mathfrak i \mathbf v_0 \cdot \mathbf r}$.
Such an initial condition can be configured by the function `qm.conf_initial_condition_Gaussian`.