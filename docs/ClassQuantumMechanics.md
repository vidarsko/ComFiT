# Class: Quantum Mechanics

Quantum mechanics is one of the most classic examples of field theories in physics. 

In this class, we simulate quantum mechanics through the evolution of the Schrödinger equation.

```python
file: comfit/models/quantum_mechanics.py 
class: QuantumMechanics
```

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

## Tutorial 1: The spreading of a free Gaussian

Set up a one-dimensional system $x \in [-10,10]$ with an initial Gaussian wave function with a width of $1$. 

??? note "Solution"
    ```python
    import comfit as cf
    import numpy as np
    import matplotlib.pyplot as plt

    qm = cf.QuantumMechanics(1,xlim=[-10,10])

    qm.conf_initial_condition_Gaussian(position=0,width=1)
    qm.plot()
    plt.show()
    ```

Calculate the probability of measuring the particle in the interval $[10,15]$.

??? note "Solution"
    ```python
    a=10
    b=15
    interval = qm.calc_region_interval(a,b)
    prob = qm.calc_integrate_field(np.abs(qm.psi)**2,interval)

    print("The probability of measiring the particle in the interval (",a,",",b,") is",prob)
    ```
    Output:
    ```bash
    The probability of measiring the particle in the interval ( 5 , 10 ) is 3.551131161818509e-07
    ```

Evolve the wave function for a time $T=4.5$, plot the result.

??? note "Solution"
    ```python
    qm.evolve_schrodinger(45)
    qm.plot()
    plt.show()
    ```

Calculate the probability of measuring the particle in the interval $[10,15]$.

??? note "Solution"
    ```bash
    The probability of measiring the particle in the interval ( 5 , 10 ) is 0.02213823670966853
    ```

Give the particle an initial velocity of $v_0=1$ and plot the initial condition. 
Evolve the particle now for $T=4.5$ and calculate the probability of measuring the particle in the aforementionned interval.

??? note "Solution"
    ```python
    qm = cf.QuantumMechanics(1,xlim=[-10,10])
    qm.conf_initial_condition_Gaussian(position=0,width=1,initial_velocity=1)

    qm.plot()
    plt.show()

    a=5
    b=10
    interval = qm.calc_region_interval(a,b)
    prob = qm.calc_integrate_field(np.abs(qm.psi)**2,interval)

    print("The probability of measiring the particle in the interval (",a,",",b,") is",prob)

    qm.evolve_schrodinger(45)
    qm.plot()
    plt.show()

    a=5
    b=10
    interval = qm.calc_region_interval(a,b)
    prob = qm.calc_integrate_field(np.abs(qm.psi)**2,interval)

    print("The probability of measiring the particle in the interval (",a,",",b,") is",prob)
    ```

Make an animation of the evolution of the particle with an initial velocity.

??? note "Solution"
    ```python
    import comfit as cf
    import numpy as np
    import matplotlib.pyplot as plt

    qm = cf.QuantumMechanics(1,xlim=[-10,10])

    qm.conf_initial_condition_Gaussian(position=0,width=1,initial_velocity=1)

    ymax = np.max(np.abs(qm.psi)**2)

    for n in range(50):
        qm.plot(ylim=ymax)
        cf.tool_save_plot(n)
        qm.evolve_schrodinger(1)

    cf.tool_make_animation_gif(n)
    ```

Now generalize these ideas to two and three dimensions.

??? note "Solution"
    Not yet written.


## Tutorial 2: Gaussian wave packet with external potentials

Initiate the same Gaussian wave packet with the initial velocity as in the previous tutorial, but add an external potential `V_ext` given by $V_{\textrm{ext}}=0.05 x^2$.
Make an animation.

??? note "Solution"
    ```python
    import comfit as cf
    import numpy as np
    import matplotlib.pyplot as plt

    qm = cf.QuantumMechanics(1,xlim=[-10,10])

    qm.V_ext = 0.05*qm.x**2

    qm.conf_initial_condition_Gaussian(position=0,width=1,initial_velocity=1)

    ymax = np.max(np.abs(qm.psi)**2)

    for n in range(200):
        qm.plot(ylim=ymax)
        cf.tool_save_plot(n)
        qm.evolve_schrodinger(5)
        plt.pause(0.01)

    cf.tool_make_animation_gif(n)
    ```

Now extend the system so that it goes from $x \in [-10,100]$ and decrease `dx` to $0.1$.
Create a Gaussian potential barrier centered at $x=30$ with a top value of $1$ and a width of $5$. 
Run the simulation. What happens?

??? note "Solution"
    ```python
    import comfit as cf
    import numpy as np
    import matplotlib.pyplot as plt

    qm = cf.QuantumMechanics(1,xlim=[-10,100],dx=0.1)

    qm.V_ext = qm.calc_Gaussian(position=30,width=10,top=1)

    qm.conf_initial_condition_Gaussian(position=0,width=1,initial_velocity=1)

    ymax = np.max(np.abs(qm.psi)**2)

    for n in range(100):
        qm.plot(ylim=ymax)
        cf.tool_save_plot(n)
        qm.evolve_schrodinger(5)

    cf.tool_make_animation_gif(n)
    ```

    Quantum tunneling. 
