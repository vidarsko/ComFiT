# Class: Quantum Mechanics

Quantum mechanics is one of the most classic examples of field theories in physics.
The Schrödinger equation is a partial differential equation that describes how the quantum state of a physical system changes over time.

In this class, we simulate quantum mechanics through the evolution of the Schrödinger equation.

```python
file: comfit/models/quantum_mechanics.py 
class: QuantumMechanics
```

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


## Central methods

`evolve_schrodinger(number_of_steps)`: Evolves the wave function psi of the quantum mechanical system for a specified number of steps.

Parameters:
`number_of_steps(int)`: The total number of time evolution steps to perform.
`conf_initial_condition_Gaussian(position, width, initial_velocity)`: Configures the initial condition of the wave function to be a Gaussian.

Parameters:
`position` (tuple or list): Initial position of the Gaussian peak.
`width` (float or tuple/list of floats): Width of the Gaussian.
`initial_velocity` (tuple or list): Initial velocity of the Gaussian.

`conf_wavefunction(psi)`: Explicitly sets the wave function to a specified field.

Parameters:
`psi (array-like)`: The wave function to set as the initial condition.


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