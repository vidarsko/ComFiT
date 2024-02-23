# Class: Bose Einstein Condensate

There are two types of particles in the world: femions and bosons.
Whereas fermions can never occupy the same quantum state due to the Pauli Exclusion principle, the same is not true for bosons.
A Bose-Einstein condensate (BEC) is a state of matter consisting of ultra-cold bosons which undergo a phase transition at a low critical temperature in which most bosons occupy the ground state of the system.
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

## Model

The BEC is in the mean field regime described by the GPE
[^dalfovo1999theory] [^kevrekidis2007emergent]. This is a non-linear
Schrödinger equation which reads
$$i\hbar \partial_t\psi = \left[-\frac{\hbar^2}{2m} \nabla^2+ V_{ext} -\mu +g|\psi|^2 \right]\psi.$$
Here $\mu$ is the chemical potential, $m$ is the mass of the bosons, $g$
is an interaction parameter and $\hbar$ is the Planc constant. $\psi$ is
the wave function describing the condensate phase and $V_{ext}$ is an
external potential. The GPE can be obtained from the varaitional
principle [^kevrekidis2007emergent][^pitaevskiiBook]

$$
\mathfrak i \hbar \partial_t \psi = \frac{\delta K}{\delta \psi^*}
$$

with the Hamiltonian

$$
K = \int d \mathbf r \left[\frac{\hbar^2}{2m}|\nabla\psi|^2 +(V_{ext} -\mu)|\psi|^2 +\frac{g}{2}|\psi|^4 \right].
$$

We introduce dimensionless units for length $\xi = \hbar/\sqrt{m\mu}$,
time $\tau = \xi/c$ and energy $E=\eta$ and rescaling the wave function
to $\psi \rightarrow \sqrt{\frac{g}{\mu}}\psi$, in addition we include a
dissipative factor $\gamma$.
This results in the dGPE on dimensionless
form as
[^gardiner2003stochastic] [^rooney2012stochastic] [^bradley2012energy] [^skaugenUnifiedPerspectiveTwodimensional2018]
$$i \partial_t \psi = (1-\mathfrak i\gamma) \left[-\frac{1}{2}\nabla^2 + V_{ext} -1 +|\psi|^2 \right]\psi.$$
$$\partial_t \psi =-i (1-\mathfrak i\gamma) \left[-\frac{1}{2}\nabla^2 + V_{ext} -1 +|\psi|^2 \right]\psi.$$
$$\partial_t \psi =-(\mathfrak i+\gamma) \left[-\frac{1}{2}\nabla^2 + V_{ext} -1 +|\psi|^2 \right]\psi.$$
$$\partial_t \psi =(\mathfrak i+\gamma) (1+\frac{1}{2}\nabla^2) \psi - (\mathfrak i + \gamma) (V_{ext} + |\psi|^2)\psi$$
In other words

$$
\omega = (\mathfrak{i}+\gamma) (1 + \frac{1}{2}\nabla^2), \quad \omega_{f} =  (\mathfrak{i}+\gamma) (1-\frac{1}{2}\mathbf{k}^2), \quad N = - (\mathfrak{i} + \gamma) (V_{ext} + |\psi|^2)\psi
$$

The evolution of the wave function is included through the function

``` {.python language="Python"}
bec.evolve_dGPE(self,number_of_steps,method='ETD2RK')
```

Where the default method is the EDT2RK method described in the next
chapter. The Hamiltoninan and its denisty in dimensionless units is
calculated by the functions

``` {.python language="Python"}
bec.calc_hamiltonian(self)
bec.calc_hamiltonian_density(self)
```

respectively.

## General form

If we include longer range interactions, we get [^steinberg2022exploring]

$$
\mathfrak{i} \partial_t \psi = (1-\mathfrak{i} \gamma)
\left[-\frac{1}{2} \nabla^2 + V_{\textrm{ext}}-1 + |\psi|^2 - g_2 \nabla^2 |\psi|^2 +  g_4 \nabla^4 |\psi|^2\right] \psi.
$$

$$
\partial_t \psi = (\mathfrak{i} + \gamma)
\left[\frac{1}{2} \nabla^2 - V_{\textrm{ext}} + 1 -  |\psi|^2 + g_2 \nabla^2 |\psi|^2- g_4 \nabla^4 |\psi|^2 \right] \psi
$$

Splitting into linear and non-linear

$$
\omega = (\mathfrak i+\gamma) \frac{1}{2} (1+\nabla^2) \quad \omega_f=   (\mathfrak i+\gamma) (1-\mathbf{k}^2)
$$

$$
N = (\mathfrak{i} + \gamma) (-V_{\textrm{ext}}  - |\psi|^2 + g_2 \nabla^2 |\psi|^2-  g_4 \nabla^4 |\psi|^2)\psi,
$$

where $g_0$ in [^steinberg2022exploring] has been set to $1$.

## Approximation of Ground States

When doing a simulation it is often convenient to start in a
configuration that is close to the ground state. We can estimate this
ground state by noticing that the GPE dissipates energy when it is
evolved in imaginary time $t \rightarrow it$
[^minguzzi2004numerical] [^kevrekidis2007emergent]. Given an external
potential $V_{ext}$ we can therefore find an approximation to the ground
state by starting with a guess and then removing energy from the guessed
state by evolving the equations in imaginary time.

To get a guess of the ground state for a given potential $V_{ext}$ we
use the Thomas-Fermi approximation
[^dalfovo1999theory] [^kevrekidis2007emergent], where we assume that
$\psi$ is slowly varying so that we can neglect the Laplacian term.
Looking for stationary solutions to the dGPE we obtain the equation

$$
0 = (V_{ext} -1 +|\psi|^2 )\psi.
$$

This has two solutions, $\psi = 0$
and

$$
|\psi|^2 = 1-V_{ext}.
$$

In the case where $V_{ext} > 1$ there is
only one possibility namely $\psi = 0$. In the case of $V_{ext}  < 1$
both the stationary solutions exists so we need to evaluate their
energy. We can do this by considering the Hamiltonian

$$
K_{TF} = \int d\mathbf r \left[(V_{ext}  -1)|\psi|^2 + \frac{1}{2} |\psi|^4 \right].
$$

Inserting the anzats $\psi =0$ we see that this vanish, while if we put
in the anzats $|\psi|^2 = 1-V_{ext}$ we get something negative. We
therefore conclude that the Thomas-Fermi ground state is given as

$$
\psi = \begin{cases}
     0 & \text{if} \quad V_{ext}  > 1 \\
     \sqrt{1 -V_{ext}}  & \text{if}\quad V_{ext}  < 1
    \end{cases}
$$

This ground state can be initialised as

``` {.python language="Python"}
bec.conf_initial_condition_Thomas_Fermi(self)
```

Once this guess is initialized we can propagate the wave function in
imaginary time using the function

``` {.python language="Python"}
bec.evolve_relax_BEC(self,number_of_steps,method='ETD2RK')
```

Note that the potential has to be a constant during this evolution.

## Potentials

A lot of the interesting dynamics of a BEC comes about when it is
interacting with an external potential. This is included as the function

``` {.python language="Python"}
bec.V_ext(t) 
```

As default this is assumed to be time independent and returns the field

``` {.python language="Python"}
bec.V_0 
```

The potential can be changed by the function

```{.python language="Python"}
bec.conf_external_potential(self, V_ext, additive=False)
```

which can be used both to set it as a function or to set it as a constant potential depending on wheter `V_ext` is a function, a constant or an numpy array. If `additive =True` one add the constant `V_ext` to the allredy existing potential.
If `V_ext` is a function it need to be on the form

```{.python language="Python"}
def V(t)
     ...
     return ...
```

The evolver will then evaluate it using the `bec.time` variable which is updated on the run.
An example using a Gaussian stirring
potential is provided in the example folder.

To make life easier we have provided a couple of popular potentials.
The harmonic potential is provided in the function

``` {.python language="Python"}
bec.set_harmonic_potential(self,R_tf)
```

Here $R_{tf}$ is the Thomas-Fermi radius [^kevrekidis2007emergent], and
the harmonic potential takes the form $$V_H = \frac{r^2}{R_{tf}^2}.$$
The Gaussian potential is provided through the function

``` {.python language="Python"}
bec.gaussian_stirring_potential(self,size,strength,position)
```

Giving the potential

$$
V_g = g e^{-|\mathbf{r} - \mathbf{r}_p|^2/\sigma^2}.
$$

Here
$g$ is the strength, $\sigma$ is the size and $\mathbf{r}_p$ is the
position.

Much of the interesting physics happens when the potential is
time-dependent. For this one can define a function as

```python
def Func():
    ...
```

and update the external potential to this function by calling

## Hydrodynamics

The dGPE can be transformed to a hydrodynamic description of the BEC.
The first step in doing this is to introduce the Madelung transformation
$\psi = \sqrt{\rho} e^{i\theta}$, where $\rho = |\psi|^2$ is
the superfluid density [^kevrekidis2007emergent]. For $\gamma = 0$ this
density is conserved and satisfy the conservation equation

$$
\partial_t \rho + \nabla\cdot \mathbf{J}_s = 0,
$$

with the superfluid current given as

$$
\mathbf{J}_s = \Im(\psi^* \nabla \psi) = \rho \nabla \theta = \rho \mathbf{v}_s.
$$

Here the superfluid velocity is introduced as
$\mathbf{v}_s = \nabla \theta$. The superfluid current can be found
with the function

``` {.python language="python"}
bec.calc_superfluid_current(self)   
```

We can also put the Madelung transformation into the Hamiltonian to get
[^bradley2012energy] [^nore1997kolmogorov]

$$
K = \int d \mathbf r \left[\frac{1}{2}\rho v_s^2 +\frac{1}{8} \frac{|\nabla \rho|^2}{\rho} + (V_{ext}-1)\rho +\frac{1}{2}\rho^4 \right].
$$

The first term here is the kinetic energy of the condensate. To
calculate this it is convenient to introduce the density weighted
velocity $\mathbf{u} = \sqrt{\rho}\mathbf{v}_s$
[^bradley2012energy]. This have the advantage of not being singular at
the centre of the topological defects. Using this we can write the
kinetic energy as $$E_k =\int d \mathbf r \frac{1}{2}u^2.$$ The density
weighted velocity and the kinetic energy can be calculated by the
functions

``` {.python language="python"}
bec.calc_velocity(self)
bec.calc_kinetic_energy(self)
```

Further if we insert the Madelung transformation into the dGPE and do
some work we can map it into the Navier-Stockes equations
[^kevrekidis2007emergent] [^bradley2012energy]

$$
\begin{aligned}
    \partial_t \rho + \nabla\cdot(\rho \mathbf v) =2\gamma \rho (1-V_{eff}),
    \\
    \partial_t \mathbf v-\mathbf v\cdot\nabla \mathbf v =  - \nabla(V_{ext} +\rho) +\frac{\gamma}{2}\nabla^2 \mathbf v.
\end{aligned}
$$

Notice that the condensate density is only conserved
when $\gamma = 0$.

## Forces on external potential

To study how the BEC are interacting with impurities one can model the
impurity as a Gaussian potential and measure the forces that are acting
on it
[^ronning2020classical] [^astrakharchik2004motion] [^pinsker2017gaussian].
From the Erhenfest theorem the forces on the condensate from the stirrer
is given as

$$
\mathbf{F} = -\langle \nabla V_{ext}\rangle.
$$

Which
means that the force acting on the stirrer is
$\mathbf{F} = \langle \nabla V_{ext}\rangle$. Written explicitly
this is

$$
\mathbf{F} = \int d\mathbf r |\psi^2| \nabla V_{ext} = -\int d\mathbf r V_{ext}\nabla|\psi^2|
.
$$

This is calculated by the function

``` {.python language="python"}
bec.calc_force_on_external_potential(self)
```

Note that this calculates the average force on all of the external
potential.

## Quantized vortices

The topological defects in the BEC takes the form of quantized vortices.
This is because the velocity is given by the gradient of the phase
$\mathbf{v}_s = \nabla \theta$. This has the consequence that the
circulation
$$\Gamma = \oint d\mathbf{l} \cdot \mathbf{v} = \oint d\mathbf{l} \cdot \nabla \theta = \int d\theta$$
is quantized in units of $2\pi$, since the field $\psi$ is single
valued. In two dimensions this vortex is a point defect, while in three
dimensions it is a line defect.
For the BEC this method is given
by the function

``` {.python language="Python"}
bec.calc_vortex_nodes( dt_psi=None)
```

wich finds the vortex nodes in two or three dimensions. If
$\partial_t \psi$ is put in as the field `dt_psi`, then the velocity of
the defects are also found.

The defects can be created in multiple ways. In addition to putting them
in by hand, one can create them by
stirring the condensate with a potential or by relaxing a disordered
initial state.

## Comoving frame

When studying a BEC that is stirred by a potential it is in some cases
convinient to consider the potential as stationary with the BEC flowing
around. This can be done by boosting the dGPE so that it reads

$$
\partial_t \psi = \mathbf V_p \cdot \nabla \psi +(\mathfrak i+\gamma) (1+\frac{1}{2}\nabla^2) \psi - (\mathfrak i + \gamma) (\mathcal U + |\psi|^2)\psi,
$$

where $\mathbf{V_p}$ is the velocity of the boost. Note that a Gallilean
boost of the GPE is often accompanied by a phase shift of the wave
function
$\psi \rightarrow \psi \exp{(\mathfrak i\mathbf V_p \cdot \mathbf r + \frac i 2 V_p^2 t)}$
which transforms the superfluid velocity to the new reference frame
[^Pismen], leaving the GPE unchanged after the Gallilean transformation.
However the equation with $\gamma \neq 0$ is not Gallilean invariant and
we therefore do not include the phase factor in the transformation. This
have the consequence that the superfluid velocity obtained from $\psi$
is measured in the lab frame and not the comoving frame
[^ronning2020classical].

To reduce the recycling of excitations into the incoming flow we
introduce a buffer region around the computational domain where $\gamma$
is large, similar to the one used in Ref. [^reeves2015identifying]. The
dissipative factor becomes a function of space and is given by
$\gamma(\mathbf{r}) =
\max[\gamma_x(x),\gamma_y(y),\gamma_z(z)]$ in three dimensions and
$\gamma(\mathbf{r}) =
\max[\gamma_x(x),\gamma_y(y)$ in two dimensions. Here

$$
\begin{aligned}
\gamma_x(x)= \frac{1}{2}\big(2 + \tanh{[(x-w_x)/d]}
-\tanh{[(x+w_x)/d]}\big) + \gamma_0.
\end{aligned}
$$

and similar for $\gamma_y(y)$ and $\gamma_z(z)$. The
constant $\gamma_0$ is the value in the bulk, $d$ is the size of the
interface between the bulk and the buffer and $w_x$ is the distance from
the centre of the domain to the buffer in the $x$ direction. Note that
when $\gamma$ is a function of space we can no longer put it into the
linear differential operator $\omega$, and we have to move it into the
non-linear part. This is taken care of in the evolver

```python
bec.evolve_comoving_dGPE(self, number_of_steps, velx,method='ETD2RK')
```

Here it is assumed that the boost is in the $x$-direction, and that the
dissipative factor is spatially dependent.


## Tutorial 1: BEC Basic framework
In this notebook we are going to discuss some basic aspects of the framework, using the Bose-Einstein condensate as an example. We start by initialising a system.



```python
import comfit as cf
import matplotlib.pyplot as plt
import numpy as np

bec = cf.BoseEinsteinCondensate(2,xRes=101,yRes=101,gamma=0.01,dt=0.1)

### setting the system to a homogeneous condensate. This will be discussed more in the next tutorial 
bec.conf_initial_condition_Thomas_Fermi()

```

Here gamma is a dissipative factor that is removing energy from the system. The last line initialised the system as a homogeneous condensate, defining the two fields bec.psi and bec.psi_f wich is the wave function in real and fourier space. 

Since the system have periodic boundary conditions we can calculate derivatives in fourier space. We have

$$ \partial_x f(r) -> i k_x f_f(k). $$ 

In the module the wave vectors k is given as the variable bec.k, and the x-component is bec.k[0]. The combination 1j * bec.k is provided in the bec.dif list. You can therefore find the x derivative of the field f in fourier space by running bec.dif[0]*f_f 




```python
#### task 1: Find the paritial derivative of psi wrt x by differentiating in fourier space
#    (hint: use np.fft.ifft2() and bec.psi_f)

```

The module contains functions for plotting different types of fields (see previous notebook), like plot_field(...) and plot_complex_field(...). 

Here is an example of how to cal a plotting function.
bec.plot_field(np.abs(bec.psi)**2,cmap_symmetric=False,colormap = 'winter')
Remember to set cmap_symmetric = False when plotting
something that is not symmetric around 0.


```python
#### task 2: Plott the derivative you found above. It will look a bit strange.

```

We also need to evolve the wave function in time. This is done by the evolver functions. The bec has three different evolvers.

bec.evolve_relax( number_of_steps, method='ETD2RK') relaxes the system by integrating the Gross-Pitaevskii equation in imaginary time. 

bec.evolve_dGPE( number_of_steps, method='ETD2RK') evolves using the damped Gross-Pitaevskii equation and

bec.evolve_comoving_dGPE(number_of_steps, velx, method='ETD2RK') evolves the system in the frame moving at speed velx in the x direction (relative to the labframe). This solver allowes for gamma to be spatialy dependent.

All the evolvers allowes you to choose which solver you want to use. The default is ETD2RK which is a second order solver. The other implemented choise is ETD4RK which is fourt order. For details see the documentation or bully Vidar. 




```python
#### task 3: Evolve the system for 50 timesteps and plot the absolut value of psi.


```

We now want to look at how to initialise and track topological defects (vortices) in the bec. We can insert a vortex dipole by running the function conf_insert_vortex_dipole(self, dipole_vector=None, dipole_position=None). If the parameters dipole_vector and dipole_position is not given the dipole is put in the midle of the domain. After this one should relax the system with evolve_relax(...) to get the right core structure. 




```python
### task 4: initialise a dipole at your favoritt position and relax the system for 10 timesteps. Plot the absolute 
# value of psi


```

Tracking of defects can be done using the function calc_vortex_nodes(self, dt_psi=None) that finds the position and charge of the defects. If dt_psi is provided it will also find the defects' velocity. The output of this function is an array of vortex dictionaries.

In the vortex dictionary the vortex's position is saved under the key word 'position', the charge under 'charge' and the velocity under 'velocity'. 

Once this array is found the vortices can be ploted using the function 
plot_vortex_nodes(self, vortex_nodes, ax=None), where vortex_nodes is the array of dictionaries that where found with calc_vortex_nodes(self, dt_psi=None) 




```python
### task 5: evolve the system 50 time-steps.

### task 6: find dt_psi (hint use the evolver a few time steps and the variable bec.dt).
# Use this to find the defects and their velocity. Using the function bec.calc_vortex_nodes( dt_psi) 
#Plot the defect positions on top of a plot of the absolute value squared of the wavefunction
#(the colormap 'gray' look very nice for this)


### task 7: try making gamma larger and see what effect this has 
```

When simulating the BEC there is a range of different properties that might be of interest. Some of them are included in the library. Here is a list of the ones that are currently included

bec.calc_superfluid_current() - returns the superfluid current

bec.calc_velocity() - returns the weighted superfluid velocity ($\vec v_s\sqrt{\rho}$) 

bec.calc_hamiltonian_density() - returns the hamiltonian density

bec.calc_hamiltonian() - returns the total energy

bec.calc_kinetic_energy() - returns the total kinetic energy ($\frac{1}{2}\int d\vec r \rho v_s^2$)

bec.calc_force_on_external_potential()- returns the total force that the condensate is exerting on the external potential. Relevant if you want to find the drag on a stirring potential. 



```python
### task 8. find and plot the superfluid current and the hamiltonian density.  
```

??? note "Solution"
    ```python
    import sys
    sys.path.append('/Users/jonasronning/Documents/Work/Numerics/Old Projects/ComFiT')


    import comfit as cf
    import matplotlib.pyplot as plt
    import numpy as np

    bec = cf.BoseEinsteinCondensate(2,xRes=101,yRes=101,gamma=0.01,dt=0.1)

    ### setting the system to a homogeneous condensate 
    bec.conf_initial_condition_Thomas_Fermi()

    #### task 1 and 2
    dphidx = np.fft.ifft2(bec.dif[0]*bec.psi_f)

    bec.plot_complex_field(dphidx)
    plt.show()
    ```


        
    ![png](Solutions%202.0_files/Solutions%202.0_0_0.png)
        



    ```python
    # task 3
    bec.evolve_dGPE(50)

    bec.plot_field(np.abs(bec.psi)**2,cmap_symmetric=False,colormap = 'winter')
    plt.show()

    ```

        evolving the dGPE: 100%|███████████████████████| 50/50 [00:00<00:00, 693.90it/s]
        


        
    ![png](Solutions%202.0_files/Solutions%202.0_1_1.png)
        



    ```python
    #Task 4
    bec.conf_insert_vortex_dipole( )

    bec.evolve_relax( 10) 

    bec.plot_field(np.abs(bec.psi)**2,cmap_symmetric=False,colormap = 'winter' )
    plt.show()
    ```

        Relaxing the BoseEinsteinCondensate...
        

        evolving the dGPE: 100%|███████████████████████| 10/10 [00:00<00:00, 398.32it/s]
        


        
    ![png](Solutions%202.0_files/Solutions%202.0_2_2.png)
        



    ```python
    # task 5 and 6

    bec.evolve_dGPE(100)
    psi_old = bec.psi
    N=10
    bec.evolve_dGPE(N)
    dt_psi = (bec.psi-psi_old)/(N*bec.dt)



    nodes = bec.calc_vortex_nodes(dt_psi) 

    ax=bec.plot_field(np.abs(bec.psi)**2,cmap_symmetric=False,colormap = 'gray')
    bec.plot_vortex_nodes(nodes,ax)
    plt.show()
    ```

        evolving the dGPE: 100%|█████████████████████| 100/100 [00:00<00:00, 662.31it/s]
        evolving the dGPE: 100%|███████████████████████| 10/10 [00:00<00:00, 670.51it/s]
        

        (101, 101)
        (101, 101)
        (101, 101)
        


        
    ![png](Solutions%202.0_files/Solutions%202.0_3_2.png)
        



    ```python
    # task 8

    J_s = bec.calc_superfluid_current()
    H = bec.calc_hamiltonian_density()

    ax=bec.plot_field(H,cmap_symmetric=False,colormap = 'cool')
    bec.plot_vector_field(J_s,ax)
    plt.show()
    ```

        0.4617129485275939
        


        
    ![png](Solutions%202.0_files/Solutions%202.0_4_1.png)
        



    ```python

    ```


## Tutorial 2: comoving frame and defect tracking

Here we are going to show how the BoseEinsteinCondensate can be evolved in the comoving frame. We start by initialising a BoseEinsteinCondensate with a gaussian potential a bit forward in the computational domain. The ground state of this configuration is found by initialising the Thomas-Fermi ground state and evolving this in imaginary time.  

A constant (in time) potential is implemented in the code as the field bec.V0. The default value is 0, i.e a homogeneous condensate. We provide two functions that calculates some commonly used potentials. The gaussian and the harmonic potential. They are given by the functions:

bec.calc_gaussian_stirring_potential(size, strength, position) and
bec.calc_harmonic_potential(R_tf)

where R_tf is the Thomas-Fermi radius (the size of the condensate, see documentation for details). 


```python
import numpy as np
import matplotlib.pyplot as plt
import comfit as cf

bec = cf.BoseEinsteinCondensate(2,xRes=256,yRes=128,gamma=0,dt=0.1)

### task 1: Set the potiential to a constant gaussian at this position [bec.xmid+50,bec.ymid] with size = 5 and
# strength = 4
# sjekk documentation



```

We now need to initialise the wave function. A convinient starting point is the ground state of the system with the given potential. Since we don't have any analytic expression for this state we find it by starting out from a guess and remove energy by evolving in imaginary time (see previous notebook and documentation). The first guess is usualy the Thomas-Fermi ground state which is implemented as the function:

bec.conf_initial_condition_Thomas_Fermi()



```python
### task 2: Initialise the wavefunction using the Thomas-Fermi ground state and relax the system in imaginary time 
# for 50 time steps. Plot the absolut value squared of the wave function
#(Hint: the evolvers vere discussed in the previous notebook)


bec.plot_field(np.abs(bec.psi)**2,cmap_symmetric=False,colormap = 'winter')
plt.show()
```

Now we sett bec.gamma to be zero inside the computational domain and 1 on the edges. This can be done with the function bec.conf_dissipative_frame(d=7, wx=50, wy=50, wz=50). Here wx/wy/wz is the distance from the midle to the start of the dissipative frame in the x/y/z direction. In two dimensions wz is not used. d is the size of the interface between the bulk and the dissipative frame. 




```python
#### task 3. make a dissipative frame with d=7, wx = 100 and wy = 50. plot bec.gamma



```

Now we want to evolve the system in the comoving frame and we want to use the function bec.calc_vortex_nodes() to find the vortices.


```python
### task 4. evolve the system in the comoving frame with vel_x = 0.4 untill t = 300. plot the absolut value squared 
# of the wavefunction and track the defects.   





```


```python
### We can now plot the vortices you traced 

ax=bec.plot_field(np.abs(bec.psi)**2,cmap_symmetric=False,colormap = 'gray')
bec.plot_vortex_nodes(nodes,ax)
plt.show()
```


```python

```




```python

```


```python

```


```python

```


??? note "Solution"
    ```python
    import numpy as np
    import matplotlib.pyplot as plt
    import comfit as cf

    bec = cf.BoseEinsteinCondensate(2,xRes=256,yRes=128,gamma=0,dt=0.1)

    ### task 1 and 2


    pot = bec.calc_gaussian_stirring_potential(5, 4, [bec.xmid+50,bec.ymid] )
    bec.conf_external_potential(pot, additive=False)

    bec.conf_initial_condition_Thomas_Fermi()

    bec.evolve_relax(50, method='ETD2RK') 

    bec.plot_field(np.abs(bec.psi)**2,cmap_symmetric=False,colormap = 'winter')
    plt.show()
    ```


    ```python
    # task 3
    bec.conf_dissipative_frame(wx=100,wy=50)

    bec.plot_field(bec.gamma,cmap_symmetric=False,colormap = 'winter')
    plt.show()



    ```


    ```python
    #task 4
    vel_x = 0.40


    t_max = 300
    timesteps = int(t_max/bec.dt)

    bec.evolve_comoving_dGPE(timesteps,vel_x,method='ETD4RK')
    N=10
    psi_old = bec.psi
    bec.evolve_comoving_dGPE(N,vel_x)
    dt_psi = (bec.psi-psi_old)/(N*bec.dt)

    nodes = bec.calc_vortex_nodes(dt_psi)
    ```


    ```python
    ax=bec.plot_field(np.abs(bec.psi)**2,cmap_symmetric=False,colormap = 'gray')
    bec.plot_vortex_nodes(nodes,ax)
    plt.show()
    ```


    ```python

    ```


## Tutorial 3: Time dependent potentials

In this notebook we are going to illustrate how the potential works in the BoseEinsteinCondensate module. We start by initialising a 2 dimensional BoseEinsteinCondensate with dimensions 101x101. The factor bec.gamma is here set to be non-zero to introduce some dissipation to the model. This is important because the bec with $\gamma = 0$ conserves energy and will start misbehaving.



```python
import comfit as cf
import matplotlib.pyplot as plt
import numpy as np

### Task 1: initialise a bec in two dimensions with resolution 101 in x and y. Make gamma = 0.05 


```

Now we need to initialize the wavefunction. Before we do that we need to specify the potential. The potential is by default a function given as

self.V_ext = lambda: self.V0

The potential can be changed by using the function bec.conf_external_potential(pot, additive=False). Here pot could be either a function or a constant



```python
### First we set the size of the harmonic
R_tf = 40

### Here we set the size and velocity of the stirrer
stirrer_radius = 20
stirrer_velocity = 0.6
freq = stirrer_velocity/stirrer_radius
size =4
strength = .9

### Defining the function for the time-dependent potential
def V_t(t):
    pos_x = bec.xmid + stirrer_radius * np.cos(freq * t)
    pos_y = bec.ymid + stirrer_radius * np.sin(freq * t)
    stirrer = bec.calc_gaussian_stirring_potential(size, strength, [pos_x, pos_y])
    harmonic = bec.calc_harmonic_potential(R_tf)
    return   harmonic + stirrer


```


```python
### Task 2: Set the potential to the time dependent one value, initialise the Thomas Fermi 
# ground state and relax the system using the  evolve_relax(...) solver for 20 time steps

# 



bec.plot_field(np.abs(bec.psi)**2,cmap_symmetric=False,colormap = 'winter')
plt.show()
```


```python




### Task 3: Evolve the system with the time-dependent potential using the ETD4RK scheme


### Task 4: Track the defects and their velocity and plot the result 


```

Now we set the potential to be time independent and run the system again. The non-zero bec.gamma are going to relax the system.


When working with a time dependent sysytem it is nice to make some movies. To do this one needs to use two functions. The first one is cf.tool_save_plot(n) wich saves the plot and label it as n. When all the plots is saved you can cal cf.tool_make_animation(N-1) which takes the figures labeled 0 - (N-1) and makes a plot. It also deletes the figures. The procedure for making a movie is therefore

for i in range(N):

    evolve(...)
    
    make_plot(...)
    
    cf.tool_save_plot(i)
   
cf.tool_make_animation(i) (notice indent)




```python
#### task 5. make an animation of the stirring potential. Evolve 10 or 20 timesteps between each frame for a total 
#### of 3000 or more timesteps.
#### Display both the absolute value squared of the wavefunction and track the vortices.  Notice that in
###  order to make the plots apear in the same axes you need to use:
###  ax = bec.plot_field(...)
###  bec.plot_vortex_nodes(nodes, ax) 
```


```python
const_pot = V_t(bec.time)


bec.conf_external_potential(const_pot, additive=False)
```


```python
timesteps = 300
bec.evolve_dGPE(timesteps,'ETD4RK')

timesteps = int(200/bec.dt)
bec.plot_field(np.abs(bec.psi)**2,cmap_symmetric=False,colormap = 'winter')
plt.show()
```

Task 6 (optional): Do the task again, but implement your own time-dependent potential.  


```python

```

??? note "Solution"
    ```python
    import comfit as cf
    import matplotlib.pyplot as plt
    import numpy as np

    ### Task 1: initialise a bec in two dimensions with resolution 101 in x and y. Make gamma = 0.05 
    bec = cf.BoseEinsteinCondensate(2,xRes=101,yRes=101,gamma=0.05,dt=0.1)

    ### First we set the size of the harmonic
    R_tf = 40

    ### Here we set the size and velocity of the stirrer
    stirrer_radius = 20
    stirrer_velocity = 0.6
    freq = stirrer_velocity/stirrer_radius
    size =4
    strength = .9

    ### Defining the function for the time-dependent potential
    def V_t(t):
        pos_x = bec.xmid + stirrer_radius * np.cos(freq * t)
        pos_y = bec.ymid + stirrer_radius * np.sin(freq * t)
        stirrer = bec.calc_gaussian_stirring_potential(size, strength, [pos_x, pos_y])
        harmonic = bec.calc_harmonic_potential(R_tf)
        return   harmonic + stirrer





    ```


    ```python
    ######## task 2 #########
    bec.conf_external_potential(V_t, additive=False)

    bec.conf_initial_condition_Thomas_Fermi()

    bec.evolve_relax(50) 

    bec.plot_field(np.abs(bec.psi)**2,cmap_symmetric=False,colormap = 'winter')
    plt.show()
    ```


    ```python
    #### task  3, 4 and 5


    bec.evolve_dGPE( 30, method='ETD4RK') 

    nodes = bec.calc_vortex_nodes()

    ax=bec.plot_field(np.abs(bec.psi)**2,cmap_symmetric=False,colormap = 'gray')
    bec.plot_vortex_nodes(nodes,ax)
    plt.show()
    ```


    ```python
    ### here is some code for making a video
    N = 300
    for n in range(N):
        bec.evolve_dGPE(10)
        ax=bec.plot_field(abs(bec.psi)**2,colormap='gray',cmap_symmetric=False,
                    clims=[0,1])
        nodes = bec.calc_vortex_nodes()
        bec.plot_vortex_nodes(nodes,ax)
        cf.tool_save_plot(n)
        

    cf.tool_make_animation(n)
    ```


    ```python
    const_pot = V_t(bec.time)


    bec.conf_external_potential(const_pot, additive=False)

    timesteps = int(200/bec.dt)
    bec.evolve_dGPE(timesteps,'ETD4RK')


    bec.plot_field(np.abs(bec.psi)**2,cmap_symmetric=False,colormap = 'winter')
    plt.show()
    ```


    ```python

    ```



## Tutorial 4: 3D BoseEinsteinCondensate in a comoving frame

We now want to take the step into 3D. This is not to different from the 2D systems.  We start by setting the potential to a gaussian and initialise the wave function by relaxing the Thomas-Fermi ground state.

```python

import numpy as np
import matplotlib.pyplot as plt
import comfit as cf

### task 1. Initialise a BoseEinsteinCondensate in 3 dimensions with resolution 64x64x64. Sett gamma = 0 and dt =0.1


### task 2. set the potential to a gaussian placed at the centre with size = 2 and strenght = 4. Initialise the
# wave function using the thomas-fermi ground-state and relax the system for 100 time steps. Plot the result

```

??? note "Solution"
    ```python
    import numpy as np
    import matplotlib.pyplot as plt
    import comfit as cf


    bec = cf.BoseEinsteinCondensate(3,xRes=64,yRes=64,zRes=64,gamma=0,dt=0.1)

    pot= bec.calc_gaussian_stirring_potential(2,4,[bec.xmid,bec.ymid,bec.zmid])
    
    bec.conf_external_potential(pot, additive=False)
    
    bec.conf_initial_condition_Thomas_Fermi()
    bec.evolve_relax(100)

    bec.plot_field(np.abs(bec.psi)**2)
    plt.show()

    ```


    ```python
    bec.conf_dissipative_frame(wx=25,wy=25,wz=25)

    bec.plot_field(bec.gamma)
    plt.show()
    ```


    ```python
    bec.psi += (0.01*np.random.randn(bec.xRes,bec.yRes,bec.zRes)+ 0.01*np.random.randn(bec.xRes,bec.yRes,bec.zRes)*(1j))*np.abs(bec.psi)**2
    bec.psi_f = np.fft.fftn(bec.psi)
    vel_x = 0.8

    t_max = 300

    timesteps = int(t_max/bec.dt)


    bec.evolve_comoving_dGPE(timesteps,vel_x,method='ETD4RK')


    nodes =  bec.calc_vortex_nodes()    

    #ax=bec.plot_field(np.abs(bec.psi)**2,cmap_symmetric=False)
    bec.plot_vortex_nodes(nodes)
    plt.show()
    ```


    ```python

    ```



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
[^steinberg2022exploring]: Steinberg, A. B., Maucher, F., Gurevich, S. V. and Thiele, U. (2022). Exploring bifurcations in Bose--Einstein condensates via phase field crystal models. Chaos: An Interdisciplinary Journal of Nonlinear Science. 32, 11 [https://doi.org/10.1063/5.0101401](https://doi.org/10.1063/5.0101401)
