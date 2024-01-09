# Base system 

This class simply initiates a system, defines the grid and contains the basic functionality for evolving in time. 

## General parameters

| Definition         | Code            | Description |
| ------------------ | --------------- | ----------- |
| $a_0$              | `a_0 = 1`       | A length scale associated with the system, in units of which all plots will be scaled. |
| $x_{\textrm{Res}}$ | `xRes = 101`    | The resolution in the x-direction. Similarly for $y_{\textrm{Res}}$, $z_{\textrm{Res}}$. |
| $\Delta x$         | `dx = 1`        | The discretization step in the x-direction. Similarly for $\Delta y$, $\Delta z$. |
| $x$                | `x`             | A $x_{\textrm{Res}}$ long array consisting of the positions. Similarly for $y$, $z$. |
| $x_{\textrm{mid}}$            | `xmid`          | The mid x-value. Similarly for $y_{\textrm{mid}}$, $z_{\textrm{mid}}$. |
| $x_{\textrm{max}}$            | `xmax`          | The size of the domain in the x-direction. Note that `xmax = x[-1] + dx`. Similarly for $y_{max}$,  $z_{max}$. |
| $\mathbf r_{\textrm{mid}}$ | `rmid`      | |
| $\Delta  t$        | `dt = 0.1`      | The time step. |

*Default values are shown as the `=value` in the code column.*


## Types of functions

There are five different types of functions:

1. `conf_`-functions: Configures the state of the system, for instance by setting an initial condition. Returns nothing.
2. `evolve_`-functions: Evolves the state in time according to some equation of motion.
3. `calc_`-functions: Calculates something from the state, returns that which has been calculated.
4. `plot_`-functions: Functions tailored to plot specific things. Returns the axes and figure.
5. `get_`-functions: Functions that return a component of a tensor field. Relevant for symmetric and antisymmetric tensors where it is not convenient to save all elements.

## Coarse-graining

A common and useful method is that of coarse-graining, which is defined as

$$
\rho = \langle \tilde \rho \rangle 
\equiv 
\int d^d r' \mathcal K(\mathbf r-\mathbf r') 
\tilde \rho(\mathbf r') ,
$$

where $\mathcal K(\mathbf r'-\mathbf r)$ is a Gaussian kernel given by 

$$
\mathcal K(\mathbf r- \mathbf r') = \frac{1}{(2\pi w^2)^{d/2}} \exp\left (-\frac{(\mathbf r-\mathbf r')^2}{2w^2}
\right ),
$$

From a numerical point of view, this is done in Fourier space since, by the convolution theorem, 

$$
\mathcal F \rho = \mathcal F {\mathcal K} \mathcal F {\tilde \rho}.
$$

Thus, we need the Fourier transform of $\mathcal K$, which is 

$$
\mathcal F {\mathcal K} = \int d^d r e^{-i \mathbf k \cdot \mathbf r} \frac{1}{(2\pi w^2)^{d/2}} \exp\left (-\frac{\mathbf r^2}{2w^2} \right )
$$
$$
= \frac{1}{(2\pi w^2)^{d/2}} \prod_{n=1}^d \int dr_n e^{-\frac{1}{2 w^2} r_n^2 -i k_n r_n} 
$$
$$
= \frac{1}{(2\pi w^2)^{d/2}} \prod_{n=1}^d \int dr_n e^{-\frac{1}{2 w^2} (r_n^2 + 2 i w^2 k_n r_n)}  
$$
$$
= \frac{1}{(2\pi w^2)^{d/2}} \prod_{n=1}^d e^{-\frac{1}{2} w^2 k_n^2} \int dr_n e^{-\frac{1}{2 w^2} (r_n + i w^2 k_n)^2} 
$$
$$ 
= e^{-\frac{1}{2} w^2 \mathbf k^2}.
$$

This is why we have the following function

```python
calc_gaussfilter_f
```
which calculates $\mathcal F [\mathcal K]$.

Typically, a field is coarse-grained with a width using the following piece of code

```python
field = np.fft.ifftn(np.fft.fftn(field) * self.calc_gaussfilter_f(width))
```

## Vortex fields

A general feature that will be used again and again is that of angle fields of vortices. 
An angle field is a field where each point in space corresponds to an angle $\theta \in \mathcal S^n$.
A vortex is a topological defect in an angle field, around which the circulation is some integer multiple of the covering of $\mathcal S^n$. 

### Angle field of a single vortex in two dimensions

In two dimensions, the angle field takes values $\theta \in [-\pi,\pi \rangle$ and a vortex is a point $\mathbf r_0$.
The angle field produced by the vortex has a circulation which is a multiple integer of $2\pi$, i.e.

$$
\oint d\theta = 2\pi s_n, 
$$

where $s_n$ is the integer charge of the vortex. 
A possible angle field for a vortex positioned at $(x_0,y_0)$ is given by 

$$
\theta_n = s_n \textrm{atan2}(y-y_0,x-x_0)
$$

### Angle field of a vortex ring in three dimensions

For a ring vortex in three dimensions centered at $\mathbf{r_0}$ with radius $R$, an angle field with the correct topological charge is given by first calculating the auxiliary quantities

$$
m_2 = (\mathbf{r}-\mathbf{r_0})\cdot \mathbf{n},
$$

$$
m_1 = |(\mathbf{r}-\mathbf{r_0}) - m_2\mathbf{n}|.
$$

and then calculating

$$
\theta_1 = \textrm{atan2}\left (m_2,m_1+R\right ) 
$$

$$
\theta_2 = \textrm{atan2}\left (m_2,m_1-R\right ) .
$$

These expressions are based on the geometry depicted in the following figure. (Note: The actual figure can't be rendered in Markdown, so please refer to the original document or image file for the figure.)

The angle field is then given by

$$
\theta(\mathbf{r}) = \textrm{mod}(\theta_1+\theta_2,[-\pi,\pi \rangle)
$$

and is implemented in the function `calc_angle_field_vortex_ring`.

### Periodic boundary conditions: Numerical implementation of angle fields

Apart from the angle field of a single vortex, the other fields are compatible with periodic boundary conditions. The expressions for these fields, however, are really only valid for an infinite region. When this is imposed on periodic boundary conditions, it results in spurious boundary effects, especially if either of the vortices is placed near the edge of the simulation domain. By simply inserting the vortices directly, we get what is shown in the following figure (a). (Note: The actual figure can't be rendered in Markdown, so please refer to the original document or image file for the figure.)

This field is not periodic on the domain. This typically causes the unintentional nucleation of vortices and strain on the boundary. We therefore seek to modify the fields so that they don't "see" the periodic boundary conditions.

In order to produce a field that is periodic on the domain, we transform the field $\theta$ to a complex field $\eta = e^{i \theta}$. The argument of this complex field has the correct winding configuration of the vortex dipole. However, we want to smooth this field so that it goes to $1$ ($\theta=0)$ at the borders. To do so, we introduce the filter function

$$
F = \frac{1}{2} \left ( \tanh((r^2-R^2)/w^2) - 1 \right ), 
$$

where $r^2 = (x-x_{\textrm{mid}})^2 + (y-y_{\textrm{mid}})^2$, which is a function that is zero in the center region and goes to $1$ at infinity over a width of $w$. The filtered field $\eta$ is then obtained by making a smooth function that goes from the previous angle field according to

$$
\tilde \eta = F \cdot \eta + (F-1) \cdot 1.
$$

$\tilde \eta$ is shown in Figure (c). The value of $w$ and $R$ can be adjusted and are found in the source code as `width` and `radius`, respectively.

From this section, it is clear that the initial vortex dipole should not be too large. Thus, we have included a warning in case it is attempted to initiate a dipole with a larger distance than a certain threshold.


## Numerical integration scheme

The systems of interest for this code are those that can be written on the form

$$
\partial_t \psi = \omega \psi + N
$$

where $\omega$ is a linear differential operator and $N$ is a non-linear operator (function of $\psi$). 
The following table shows some examples from the models that we will discuss in the following chapters.

| Model | $\omega$ | $\omega_f(\mathbf{k})$ | $N$ |
| --- | --- | --- | --- |
| Quantum Mechanics | $\frac{1}{2}i \nabla^2 $ | $-\frac{1}{2} i \mathbf{k}^2$ | $- i V$ |
| BEC | $(i+\gamma) (1+\frac{1}{2}\nabla^2)$ | $(i+\gamma) (1-\frac{1}{2}\mathbf k^2)$ | $- (i + \gamma) (V_{ext} + \psi \psi^*)\psi$ |
| Active Nematic | $\frac{K}{\gamma} \nabla^2 +\frac{AB}{\gamma}$ | $-\frac{K}{\gamma} k^2 +\frac{AB}{\gamma}$ | $- \mathbf{u}\cdot \nabla Q + Q \Omega -\Omega Q - \frac{2A}{\gamma}Q^2_{kk}Q$ |

*Table: Examples of time evolution operators, non-conserved.*

In the following, we will explain the method of evolution of exponential time differencing as introduced in Ref. [coxExponentialTimeDifferencing2002](References.md) for stiff systems. This will result in two integration schemes, the exponential time differencing second order Runge Kutta 2 (ETD2RK) scheme and the forth order ETD4RK scheme. As in Ref. [coxExponentialTimeDifferencing2002], we will show an intuitive way to obtain the former and only recite the expressions for the latter.

### The ETD2RK scheme

To see how we obtain the ETD2RK scheme, we take the Fourier transformation of Eq. `ref:timeevolution`, and get:

$$
\partial_t \psi_f = \omega_f \psi_f + N_f .
$$

$$
(\partial_t \psi_f)e^{- \omega_f t}  -  \omega_f \psi_f  e^{- \omega_f t} = N_f e^{- \omega_f t} .
$$

$$
\partial_t (\psi_f e^{-\omega_f t} ) =  e^{-\omega_f t} N_f
$$

Where $\psi_f(\mathbf{k})$ is the Fourier transform of $\psi(\mathbf{r})$. Integrating from $t$ to $t+ \Delta t$, one finds:

$$
\psi_f (t+\Delta t) e^{- \omega_f(t+ \Delta t)} - \psi_f (t) e^{- \omega_f t} = \int_t^{t+\Delta t} e^{- \omega_f \tau} N_f d\tau
$$

we multiply by $e^{ \omega_f(t+ \Delta t)}$ and get:

$$
\psi_f (t+\Delta t) = \psi_f (t) e^{\omega_f \Delta t} + e^{ \omega_f (t+\Delta t)} \int_t^{t+\Delta t} e^{- \omega_f \tau} N_f d\tau
$$

This is an exact result, however, the last integral is unknown. In order to calculate the last integral here, we approximate it by $\psi_f N (t+\tau) \approx N_{f0} +  \frac{\Delta N_f}{\Delta t} \tau$ where $N_{f0} = N_f(\psi(t))$ and $\Delta N_f = N_f(t+\Delta t)-N_f(t)$. We also change the integration limits from $\tau \in [t,t+\Delta t]$ to $\tau \in [0,\Delta t]$, which gives:

$$
\psi_f (t+\Delta t) = \psi_f (t) e^{ \omega_f \Delta t} $$
$$+ e^{\omega_f \Delta t} \frac{1}{- \omega_f} [e^{- \omega_f \tau}]_0^{\Delta t} N_{f0} + e^{ \omega_f \Delta t} \frac{1}{\Delta t} [\frac{\tau e^{-\omega_f \tau}}{-\omega_f} - \frac{e^{-\omega_f \tau}}{\omega_f^2}]_0^{\Delta t} \Delta N_f
$$

To find $\psi_f (t+\Delta t)$, we would need to know the value $N_f (t+\Delta t)$ before finding the state at $\psi(t+\Delta t)$. To do this, we first find a predicted state $\psi_a$ by assuming $\Delta N_f=0$ and calculating $\psi(t)$ according to the equation above. This lets us calculate an approximate $\Delta N_f = N_{fa} - N_{f0}$ and we use this in order to evolve $\psi$. This is the ETD2RK scheme.

---
$$
\psi_{fa} = \psi_{f0} + I_{f1} \psi_{f0}
$$

$$
\psi_f (t+\Delta t) = \psi_{fa} + I_{f2} (N_{fa} - N_{f0})
$$

$$
\textrm{where}
$$

$$
I_{f0} = e^{\omega_f \Delta t}
$$

$$
I_{f1} = \frac{1}{\omega_f} (e^{ \omega_f \Delta t} - 1)
$$

$$
I_{f2} = \frac{1}{\Delta t \omega_f^2} (e^{ \omega_f \Delta t} -1  -\omega_f \Delta t)
$$

---

For numerical purposes, it is useful to calculate the small $\omega_f$ limit. We expand the exponential in its Taylor series and keep the leading order term to get:

$$
I_{f0} \approx 1
$$

$$
I_{f1} \approx   \frac{1}{\omega_f} (1 + \omega_f \Delta t - 1) = \Delta t
$$

$$
I_{f2} \approx \frac{1}{\Delta t \omega_f^2}
 \left (
1 + \omega_f \Delta t + \frac{1}{2} ( \omega_f \Delta t )^2
 -1  
 -\omega_f \Delta t 
\right ) = \frac{1}{2} \Delta t
$$

In $I_{f1}$, and $I_{f2}$ there is a division by $0$ when $\omega_f = 0$. To avoid numerical issues related to this we use the above limits when $|\omega_f|$ is smaller than a tolerance. We don't use the limit for $I_{f0}$ since it doesn't contain a division by $0$. The function `evolve_ETD2RK_loop` defined in the base system class performs an ETD2RK step. This function is called by the evolvers discussed in the model chapter if the method is defined as `method = "ETD2RK"`. This is the default solver if `method` is not set. The integrating factors for a given $\omega_f(\mathbf{k})$ can be found with the function `calc_evolution_integrating_factors_ETD2RK` where the variable `tol` gives when the factors should be replaced by their leading order Taylor expansion.
Note that all solvers defined in the  class \lstinline{BaseSystem} updates the time variable
`self.t` to allow for time-dependents in the non-linear term.

### The ETD4RK scheme

Following Ref. \cite{coxExponentialTimeDifferencing2002}, we may generalize the method to a fourth order Runge-Kutta as follows

---
$$
\begin{aligned}
\psi_{fa} &= I_{f0} \psi_{f0} +  I_{f1} N_{f0} \\
\psi_{fb} &= I_{f0} \psi_{f0} + I_{f1} N_{fa} \\
\psi_{fc} &= I_{f0} \psi_{fa} + I_{f1} (2 N_{fb} - N_{f0}) \\
\psi_f (t+\Delta t) &;= I_{f2} \psi_{f0} + I_{f3} N_{f0} + I_{f4} (N_{fa} + N_{fb}) + I_{f5} N_{fc}
\end{aligned}
$$

where

$$ 
\begin{aligned}
I_{f0} &= e^{\omega_f \Delta t/2} \\
I_{f1} &= \frac{1}{\omega_f}
( e^{ \omega_f \Delta t/2} - 1) \\
I_{f2} &= e^{\omega_f \Delta t} \\
I_{f3} &= \frac{1}{ \omega_f^3\Delta t^2} 
\left (
-4 -  \Delta t \omega_f + e^{\omega_f \Delta t}(4-3\omega_f \Delta t + \omega_f^2 \Delta t^2 ) 
\right ) \\
I_{f4} &= \frac{2}{ \omega_f^3\Delta t^2}
\left (
2 + \omega_f \Delta t + e^{\omega_f \Delta t}(-2 + \omega_f \Delta t)
\right ) \\
I_{f5} &= \frac{1}{ \omega_f^3\Delta t^2}
\left (
-4 - 3 \omega_f \Delta t -  \omega_f^2 \Delta t^2 + e^{\omega_f \Delta t}(4-\omega_f \Delta t)
\right )
\end{aligned}
$$

---

*Algorithm: The ETD4RK scheme*

In the small $\omega_f$ limit, we have
$$
I_{f0} \approx 1
$$

$$
I_{f1} \approx \frac{1}{2} \Delta t
$$

$$
I_{f2} \approx 1
$$

$$
I_{f3} \approx 
\frac{1}{ \omega_f^3\Delta t^2} \times 
\left (
-4 -  \Delta t \omega_f + (1 + \omega_f \Delta t + \frac{1}{2} (\omega_f \Delta t)^2 + \frac{1}{6} (\omega_f \Delta t)^3 )(4-3\omega_f \Delta t + \omega_f^2 \Delta t^2 )
\right ) $$

$$
= \frac{1}{ \omega_f^3\Delta t^2} 
\left (
\frac{4}{6} (\omega_f \Delta t)^3 - \frac{3}{2} (\omega_f \Delta t)^3 + (\omega_f \Delta t)^3
\right )
= \frac{1}{6} \Delta t
$$

$$
I_{f4} \approx \frac{2}{ \omega_f^3\Delta t^2}
\left (
2 + \omega_f \Delta t +(1 + \omega_f \Delta t + \frac{1}{2} (\omega_f \Delta t)^2 + \frac{1}{6} (\omega_f \Delta t)^3 )(-2 + \omega_f \Delta t)
\right ) 
$$

$$
=
\frac{2}{ \omega_f^3\Delta t^2}
\left (
\frac{1}{2} (\omega_f \Delta t)^3-\frac{2}{6}(\omega_f \Delta t)^3
\right )
= \frac{1}{3} \Delta t
$$

$$
I_{f5} = 
\frac{1}{ \omega_f^3\Delta t^2} \times 
\left (
-4 - 3 \omega_f \Delta t -  \omega_f^2 \Delta t^2 + (1 + \omega_f \Delta t + \frac{1}{2} (\omega_f \Delta t)^2 + \frac{1}{6} (\omega_f \Delta t)^3 )(4-\omega_f \Delta t)
\right )  
$$

$$
=\frac{1}{ \omega_f^3\Delta t^2} 
\left (
\frac{4}{6} (\omega_f \Delta t)^3 - \frac{1}{2} (\omega_f \Delta t)^3
\right )
= \frac{1}{6} \Delta t
$$

Similar as for the EDT2RK case $I_{f1}$, $I_{f3}$, $I_{f4}$, and $I_{f5}$ contains a division by $0$ when $\omega_f = 0$.  
We therfore replace these coeficients with their limits when $|\omega_f|$ is smaller than a tolerance. 
This has been important in order to make the the code stable for some of the systems. 
In the same way as the EDT2RK scheme this is implemented as the function
`self.evolve_ETD4RK_loop(self, integrating_factors_f, non_linear_evolution_function_f, field, field_f)`
This function is called by the evolvers discussed in the model chapter if the method is defined as ```method = "ETD4RK"```, the integrating factors are found with
`self.calc_evolution_integrating_factors_ETD4RK(self, omega_f, tol=10**(-4))`.

### The fully non-linear limit
It is both interesting and enlightening to see the fully non-linear limit of these equations, i.e. the limit in which $\omega_f =0$, $N_f = \partial_t \psi \equiv \dot{\psi}_f$ and Eqs. (\ref{eq:ETD2RK_IF0_small_omega_limit}-\ref{eq:ETD2RK_IF2_small_omega_limit}) and (\ref{eq:ETD4RK_IF0_small_omega_limit}-\ref{eq:ETD4RK_IF5_small_omega_limit}) are exact. 
For the ETD2RK scheme, we get 
$$
\psi_{af} = \psi_{0f} +  \dot{\psi}_{0_f} \Delta t 
$$
$$
\psi(t+\Delta t)_f = \psi_{0f} + \dot{\psi}_{0f} \frac{\Delta t}{2} + \dot{\psi}_{af} \frac{\Delta t}{2},
$$
which is a two-stage Runge-Kutta method called Heun's method.

The ETD4RK scheme becomes
$$
\psi_{af} = \psi_{0f} +  \dot{\psi}_{0f} \frac{\Delta t}{2} 
$$
$$
\psi_{bf} =  \psi_{0f} +  \dot{\psi}_{af} \frac{\Delta t}{2} 
$$
$$
\psi_{cf} = \psi_{af} + (2\dot{\psi}_{bf} - \dot{\psi}_{0f}) \frac{\Delta t}{2} 
$$
$$
\psi (t+\Delta t)_f = \psi_{0f} + \frac{1}{6} (\dot{\psi}_{0f} + 2 \dot{\psi}_{af} + 2 \dot{\psi}_{bf} + \dot{\psi}_{cf} ) \Delta t.
$$
Note that this is not the typical Runge-Kutta 4 method, due to the differences in calculating $\psi_{cf}$.
The reason is that a straight-forward generalization of the Runge-Kutta 4 method will not produce a fourth-order method in the general case [coxExponentialTimeDifferencing2002](References.md).

### The fully linear limit

If $N=0$, the evolution equation changes to 
$$
\psi_f(t+\Delta t) = e^{\omega_f \Delta t} \psi_f.
$$
An example of this is the schrödinger equation, for which $\omega_f = -\frac{1}{2}  i \mathbf k^2$, so we get
$$
\psi_f(t+\Delta t) = e^{- i \frac{1}{2} \mathbf k^2 \Delta t} \psi_f.
$$
This is an exact eqution, of course, so you may evolve this free particle solution to any time.

# Algorithms for tracking defects 

To be written

## Calculating the velocity

The equations for the velocity are taken from Ref.
[@skogvollUnifiedFieldTheory2023](References.md), simplified using Mathematica and then
substituted for python code using chatGPT.

# Plotting

The package comes with a lot of plotting functions, so

The default plot command is simply `plot()`, which will plot the current
state of the system according to some arbitrarily chosen standard.

## Plotting with matplotlib

It is useful to recap how plotting works with python briefly. Matplotlib
operates with three levels of plotting: *figures*, *axes* and
*everything else*. The figure represents the plotting window that shows
your plot, whereas axes are the individual plots represented on the
figure. In the axes, you may place different things, like image
instances, line2D objects or text. 

Therefore, the implemented plot functions will all take the ax as an
optional input and give that as an optional output if not provided.

Each class has a function called `plot`, which is our best idea of how
to plot the current configuration of the field.

## Predefined figure sizes

Obviously, a lot of these plots are meant for publication. Therefore,
there is a tool package that contains predefined and fitting sizes for
the figures in question. $$TO BE ADDED$$

## The plot in plane function

(Essentially interpolation - need to look at the packag made for that and see if it is simpler)

## Animation




