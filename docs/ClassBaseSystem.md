# Class: Base system

This class simply initiates a system, defines the grid and contains the basic functionality for evolving in time.

## General keywords and parameters

The only required input argument to the BaseSystem class is the `dim` argument, which specifies the dimension of the system.
In some cases, default values of other parameter depend on the value of `dim`, and are represented by curly brackets:

$$
\left \lbrace \begin{array}{l} \textrm{default value if } \texttt{dim }= 1 \\ \textrm{default value if } \texttt{dim }= 2  \\ \textrm{default value if } \texttt{dim }= 3  \\ \end{array} \right \rbrace
$$

These are the optional keywords for the `BaseSystem` class.

| Keyword | Definition | Default value|
|---------|------------|--------------|
| `xmin`  | Minimum value of $x$ of the simulation domain | $0$ |
| `ymin`  | Minimum value of $y$ of the simulation domain | $0$ |
| `zmin`  | Minimum value of $z$ of the simulation domain | $0$ |
| `xmax`  | Maximum value of $x$ of the simulation domain. | $100$ |
| `ymax`  | Maximum value of $y$ of the simulation domain | $\left \lbrace \begin{array}{c} 1 \\ 100 \\ 100 \\ \end{array} \right \rbrace$ |
| `zmax`  | Maximum value of $z$ of the simulation domain | $\left \lbrace \begin{array}{c} 1 \\ 1 \\  100 \\ \end{array} \right \rbrace$ |
| `xRes`  | Resolution of the $x$ axis | $101$ |
| `yRes`  | Resolution of the $y$ axis | $\left \lbrace \begin{array}{c} 1 \\ 101 \\  101 \\ \end{array} \right \rbrace$ |
| `zRes`  | Resolution of the $z$ axis | $\left \lbrace \begin{array}{c} 1 \\ 1 \\  101 \\ \end{array} \right \rbrace$ |
| `dx`    | Spacing between points on the $x$-axis. Trumps `xRes` if provided. `xmax` will be modified to match. | $\frac{\texttt{xmax}-\texttt{xmin}}{\texttt{xRes}} = 1$ |
| `dy`    | Spacing between points on the $y$-axis. Trumps `yRes` if provided. | $\frac{\texttt{ymax}-\texttt{ymin}}{\texttt{yRes}} = 1$ |
| `dz`    | Spacing between points on the $x$-axis. Trumps `zRes` if provided. | $\frac{\texttt{zmax}-\texttt{zmin}}{\texttt{zRes}} = 1$ |
| `xlim`  | List or tuple consisting of the lower and upper limit for the simulation domain in the $x$-direction. Trumps `xmin` and `xmax` if provided. | $(\texttt{xmin},\texttt{xmax}) = (0,101)$ |
| `ylim`  | List or tuple consisting of the lower and upper limit for the simulation domain in the $y$-direction. Trumps `ymin` and `ymax` if provided. | $(\texttt{ymin},\texttt{ymax}) = \left \lbrace \begin{array}{c} (0,1) \\ (0,101) \\  (0,101) \\ \end{array} \right \rbrace$|
| `zlim`  | List or tuple consisting of the lower and upper limit for the simulation domain in the $z$-direction. Trumps `zmin` and `zmax` if provided. | $(\texttt{xmin},\texttt{xmax}) = \left \lbrace \begin{array}{c} (0,1) \\ (0,1) \\  (0,101) \\ \end{array} \right \rbrace$|
| `time` | Float specifying the time of initialization | $0$ |
| `a0` | Characteristic length scale associated with the system, in units of which all plots will be scaled. This is also the default width used with the coarse-graining operation. | $1$ |

From these keywords, a number of useful parameters are constructed, given in the table below.

| Parameter      | Definition   |
| -------------- | --------------|
| `x`            | Numpy array with dimensions $\left \lbrace \begin{array}{l} \texttt{xRes} \\ \texttt{xRes}\times 1  \\ \texttt{xRes}\times 1 \times 1  \\ \end{array} \right \rbrace$ consisting of the grid points from `xmin` to (including) `xmax-dx`. |
| `y`            | Numpy array with dimensions $\left \lbrace \begin{array}{l} 1 \\ 1 \times \texttt{yRes}  \\ 1 \times \texttt{yRes} \times 1  \\ \end{array} \right \rbrace$ consisting of the grid points from `ymin` to (including) `ymax-dy`. |
| `z`            | Numpy array with dimensions $\left \lbrace \begin{array}{l} 1 \\ 1  \\ 1 \times 1 \times \texttt{zRes} \\ \end{array} \right \rbrace$ consisting of the grid points from `zmin` to (including) `zmax-dz`. |
| `xmidi`        | Index of the mid $x$-value. In the case of an odd `xRes`, this midpoint index will not hit the middle exactly but undershoot by `dx/2`. |
| `xmid`         | The $x$ value given by `xmidi`. |
| `ymidi`        | Index of the mid $y$-value. In the case of an odd `yRes`, this midpoint index will not hit the middle exactly but undershoot by `dy/2`. |
| `ymid`         | The $y$ value given by `ymidi`. |
| `zmidi`        | Index of the mid $z$-value. In the case of an odd `zRes`, this midpoint index will not hit the middle exactly but undershoot by `dz/2`. |
| `zmid`         | The $z$ value given by `zmidi`. |

Note that even though variables like `yRes`, `zRes` etc. are defined in cases where they are not relevant, such as for a $1$-dimensional system, they play no significant role in any calculations in such situations.

![](images/base_system_x_axis_illustration.png#only-light)
![](images/base_system_x_axis_illustration-colorinverted.png#only-dark)


Periodic boundary conditions means that `xmax` and `xmin` are identified as the same point. 

## Types of functions

There are five different types of functions:

1. `conf_`-functions: Configures the state of the system, for instance by setting an initial condition. Output nothing.
2. `evolve_`-functions: Evolves the state in time according to some equation of motion.
3. `calc_`-functions: Calculates something from the state, returns that which has been calculated.
4. `plot_`-functions: Functions tailored to plot specific things. Output the axes and figure.
5. `get_`-functions: Functions that return a component of a tensor field. Relevant for symmetric and antisymmetric tensors where it is not convenient to save all elements.

## Coarse-graining

A common and useful method is that of coarse-graining, which is defined as

$$
\rho = \langle \tilde \rho \rangle
\equiv \int d^d r' \mathcal K(\mathbf r-\mathbf r')
\tilde \rho(\mathbf r') ,
$$

where $\mathcal K(\mathbf r'-\mathbf r)$ is a Gaussian kernel given by

$$
\mathcal K(\mathbf r- \mathbf r') = \frac{1}{(2\pi w^2)^{d/2}} \exp\left (-\frac{(\mathbf r-\mathbf r')^2}{2w^2}
\right ),
$$

From a numerical point of view, this is done in Fourier space since, by the convolution theorem,

$$
\rho_{\mathfrak f} = \mathcal K_{\mathfrak f} \tilde \rho_{\mathfrak f}.
$$

Thus, we need the Fourier transform of $\mathcal K$, which is

$$
\mathcal K_{\mathfrak f} = \int d^d r e^{-i \mathbf k \cdot \mathbf r} \frac{1}{(2\pi w^2)^{d/2}} \exp\left (-\frac{\mathbf r^2}{2w^2} \right )
$$

$$
= \frac{1}{(2\pi w^2)^{d/2}} \prod_{n=1}^d \int dr_n e^{-\frac{1}{2 w^2} r_n^2 - \mathfrak i k_n r_n}
$$

$$
= \frac{1}{(2\pi w^2)^{d/2}} \prod_{n=1}^d \int dr_n e^{-\frac{1}{2 w^2} (r_n^2 + 2 \mathfrak i w^2 k_n r_n)}  
$$

$$
= \frac{1}{(2\pi w^2)^{d/2}} \prod_{n=1}^d e^{-\frac{1}{2} w^2 k_n^2} \int dr_n e^{-\frac{1}{2 w^2} (r_n + \mathfrak i w^2 k_n)^2}
$$

$$
= e^{-\frac{1}{2} w^2 \mathbf k^2}.
$$

This is why we have the following function

```python
calc_Gaussian_filter_f
```

which calculates $\mathcal K_{\mathfrak f}$.

Typically, a field is coarse-grained with a width using the following piece of code

```python
field = sp.fft.ifftn(sp.fft.fftn(field) * self.calc_Gaussian_filter_f(width))
```

The Gaussian function is actually so useful that is given by can be calculated using

```python
Gaussian = bs.calc_Gaussian()
```

## Vortex fields

A general feature that will be reused is that of vortex fields.
An angle field is a field where each point in space corresponds to an angle $\theta \in \mathcal S^n$.
A vortex is a topological defect in an angle field, around which the circulation is some integer multiple of the covering of $\mathcal S^n$.

### Angle field of a single vortex in two dimensions

In two dimensions, the angle field takes values $\theta \in [-\pi,\pi \rangle$ and a vortex is a point $\mathbf r_0$.
The angle field produced by the vortex has a circulation which is a multiple integer of $2\pi$, i.e.,

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

These expressions are based on the geometry depicted in the following figure.

![Vortex ring angle field explanation](images/base_system_vortex_ring_angle_field_explanation.png#only-light)
![Vortex ring angle field explanation](images/base_system_vortex_ring_angle_field_explanation-colorinverted.png#only-dark)


*Vortex ring angle field explanation:* Geometry of a vortex ring in the plane given by $\vec n$.
$\mathcal N'$ is the plane normal to the tangent vector $\vec t'$ at $\vec r'$ upon which we impose a Cartesian coordinate system to determine the angles $\theta_1$, $\theta_2$ that are used to construct the (inset) initial angle field.
Figure reprinted from Ref.[^skogvollPhaseFieldCrystal2022] with permission.

The angle field is then given by

$$
\theta(\mathbf{r}) = \textrm{mod}(\theta_1+\theta_2,[-\pi,\pi \rangle)
$$

and is implemented in the function `calc_angle_field_vortex_ring`.

### Periodic boundary conditions: Numerical implementation of angle fields

Apart from the angle field of a single vortex, the other fields are compatible with periodic boundary conditions. The expressions for these fields, however, are really only valid for an infinite region. When this is imposed on periodic boundary conditions, it results in spurious boundary effects, especially if either of the vortices is placed near the edge of the simulation domain. By simply inserting the vortices directly, we get what is shown in the following figure (a).

![Numerical implementaiton of periodic angle fields](images/base_system_numerical_implementation_of_periodic_angle_fields.png#only-light)
![Numerical implementaiton of periodic angle fields](images/base_system_numerical_implementation_of_periodic_angle_fields-colorinverted.png#only-dark)


*Numerical implementaiton of periodic angle fields:*
The angle field of panel (a) has been filtered by the field $F$ with $w=0.2x_{\textrm{max}}$ to produce the periodic field given in panel (c).
This field is simply rolled to produce a different position for the dipole in panel (d).

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

| Model | $\omega$ | $\omega_{\mathfrak f}(\mathbf{k})$ | $N$ |
| --- | --- | --- | --- |
| Quantum Mechanics | $\frac{1}{2}i \nabla^2 $ | $-\frac{1}{2} i \mathbf{k}^2$ | $- i V$ |
| BEC | $(i+\gamma) (1+\frac{1}{2}\nabla^2)$ | $(i+\gamma) (1-\frac{1}{2}\mathbf k^2)$ | $- (i + \gamma) (V_{ext} + \psi \psi^*)\psi$ |
| Active Nematic | $\frac{K}{\gamma} \nabla^2 +\frac{AB}{\gamma}$ | $-\frac{K}{\gamma} k^2 +\frac{AB}{\gamma}$ | $- \mathbf{u}\cdot \nabla Q + Q \Omega -\Omega Q - \frac{2A}{\gamma}Q^2_{kk}Q$ |

*Table: Examples of time evolution operators, non-conserved.*

In the following, we will explain the method of evolution of exponential time differencing for stiff systems[^coxExponentialTimeDifferencing2002].
This will result in two integration schemes, the exponential time differencing second order Runge Kutta 2 (ETD2RK) scheme and the forth order ETD4RK scheme. As in Ref. [coxExponentialTimeDifferencing2002], we will show an intuitive way to obtain the former and only recite the expressions for the latter.

### The ETD2RK scheme

To see how we obtain the ETD2RK scheme, we take the Fourier transformation of Eq. `ref:timeevolution`, and get:

$$
\partial_t \psi_{\mathfrak f} = \omega_{\mathfrak f} \psi_{\mathfrak f} + N_{\mathfrak f} .
$$

$$
(\partial_t \psi_{\mathfrak f})e^{- \omega_{\mathfrak f} t}  -  \omega_{\mathfrak f} \psi_{\mathfrak f}  e^{- \omega_{\mathfrak f} t} = N_{\mathfrak f} e^{- \omega_{\mathfrak f} t} .
$$

$$
\partial_t (\psi_{\mathfrak f} e^{-\omega_{\mathfrak f} t} ) =  e^{-\omega_{\mathfrak f} t} N_{\mathfrak f}
$$

Where $\psi_{\mathfrak f}(\mathbf{k})$ is the Fourier transform of $\psi(\mathbf{r})$. Integrating from $t$ to $t+ \Delta t$, one finds:

$$
\psi_{\mathfrak f} (t+\Delta t) e^{- \omega_{\mathfrak f}(t+ \Delta t)} - \psi_{\mathfrak f} (t) e^{- \omega_{\mathfrak f} t} = \int_t^{t+\Delta t} e^{- \omega_{\mathfrak f} \tau} N_{\mathfrak f} d\tau
$$

we multiply by $e^{ \omega_{\mathfrak f}(t+ \Delta t)}$ and get:

$$
\psi_{\mathfrak f} (t+\Delta t) = \psi_{\mathfrak f} (t) e^{\omega_{\mathfrak f} \Delta t} + e^{ \omega_{\mathfrak f} (t+\Delta t)} \int_t^{t+\Delta t} e^{- \omega_{\mathfrak f} \tau} N_{\mathfrak f} d\tau
$$

This is an exact result, however, the last integral is unknown. In order to calculate the last integral here, we approximate it by $N (t+\tau) \approx N_{\mathfrak f 0} +  \frac{\Delta N_{\mathfrak f}}{\Delta t} \tau$ where $N_{\mathfrak f 0} = (N(\psi(t))_{\mathfrak f}$ and $\Delta N_{\mathfrak f} = N_{\mathfrak f}(t+\Delta t)-N_{\mathfrak f}(t)$. We also change the integration limits from $\tau \in [t,t+\Delta t]$ to $\tau \in [0,\Delta t]$, which gives:

$$
\psi_{\mathfrak f} (t+\Delta t) = \psi_{\mathfrak f} (t) e^{ \omega_{\mathfrak f} \Delta t}
$$

$$
+ e^{\omega_{\mathfrak f} \Delta t} \frac{1}{- \omega_{\mathfrak f}} [e^{- \omega_{\mathfrak f} \tau}]_0^{\Delta t} N_{\mathfrak f 0} + e^{ \omega_{\mathfrak f} \Delta t} \frac{1}{\Delta t} [\frac{\tau e^{-\omega_{\mathfrak f} \tau}}{-\omega_{\mathfrak f}} - \frac{e^{-\omega_{\mathfrak f} \tau}}{\omega_{\mathfrak f}^2}]_0^{\Delta t} \Delta N_{\mathfrak f}
$$

To find $\psi_{\mathfrak f} (t+\Delta t)$, we would need to know the value $N_{\mathfrak f} (t+\Delta t)$ before finding the state at $\psi(t+\Delta t)$. To do this, we first find a predicted state $\psi_a$ by assuming $\Delta N_{\mathfrak f}=0$ and calculating $\psi(t)$ according to the equation above. This lets us calculate an approximate $\Delta N_{\mathfrak f} = N_{\mathfrak f a} - N_{\mathfrak f 0}$ and we use this in order to evolve $\psi$. This is the ETD2RK scheme.

---
$$
\psi_{\mathfrak f a} = I_{\mathfrak f 0} \psi_{\mathfrak f 0} + I_{\mathfrak f 1} N_{\mathfrak f 0}
$$

$$
\psi_{\mathfrak f} (t+\Delta t) = \psi_{\mathfrak f a} + I_{\mathfrak f 2} (N_{\mathfrak f a} - N_{\mathfrak f 0})
$$

$$
\textrm{where}
$$

$$
I_{\mathfrak f 0} = e^{\omega_{\mathfrak f} \Delta t}
$$

$$
I_{\mathfrak f 1} = \frac{1}{\omega_{\mathfrak f}} (e^{ \omega_{\mathfrak f} \Delta t} - 1)
$$

$$
I_{\mathfrak f 2} = \frac{1}{\omega_{\mathfrak f}^2 \Delta t} (e^{ \omega_{\mathfrak f} \Delta t} -1  -\omega_{\mathfrak f} \Delta t)
$$

$$
N_{\mathfrak f 0} = (N(\psi(t),t))_{\mathfrak f}
$$

$$
N_{\mathfrak f a} = (N(\psi_a,t+\Delta t))_{\mathfrak f}
$$

---

Note that $N_{\mathfrak f}$ is a non-linear function of the field variable $\psi$, but can also be an explicit variable of time $t$, i.e. $N_{\mathfrak f}(\psi,t)$.
Therefore, in the code, it has to be encoded as a function of these two variables `calc_nonlinear_evolution_function_f(self, psi, t)`.

For numerical purposes, it is useful to calculate the small $\omega_{\mathfrak f}$ limit. We expand the exponential in its Taylor series and keep the leading order term to get:

$$
I_{\mathfrak f 0} \approx 1
$$

$$
I_{\mathfrak f 1} \approx   \frac{1}{\omega_{\mathfrak f}} (1 + \omega_{\mathfrak f} \Delta t - 1) = \Delta t
$$

$$
I_{\mathfrak f 2} \approx \frac{1}{\omega_{\mathfrak f}^2 \Delta t}
 \left ( 1 + \omega_{\mathfrak f} \Delta t + \frac{1}{2} ( \omega_{\mathfrak f} \Delta t )^2
 -1 - \omega_{\mathfrak f} \Delta t
\right ) = \frac{1}{2} \Delta t
$$

In $I_{\mathfrak f 1}$, and $I_{\mathfrak f 2}$ there is a division by $0$ when $\omega_{\mathfrak f} = 0$. To avoid numerical issues related to this we use the above limits when $|\omega_{\mathfrak f}|$ is smaller than a tolerance. We don't use the limit for $I_{\mathfrak f 0}$ since it doesn't contain a division by $0$. The function `evolve_ETD2RK_loop` defined in the base system class performs an ETD2RK step. This function is called by the evolvers discussed in the model chapter if the method is defined as `method = "ETD2RK"`. This is the default solver if `method` is not set. The integrating factors for a given $\omega_{\mathfrak f}(\mathbf{k})$ can be found with the function `calc_evolution_integrating_factors_ETD2RK` where the variable `tol` gives when the factors should be replaced by their leading order Taylor expansion.
Note that all solvers defined in the  class \lstinline{BaseSystem} updates the time variable
`self.t` to allow for time-dependents in the non-linear term.

### The ETD4RK scheme

Following Ref.[^coxExponentialTimeDifferencing2002], we may generalize the method to a fourth order Runge-Kutta as follows

---
$$
\begin{aligned}
\psi_{\mathfrak f a} &= I_{\mathfrak f 0} \psi_{\mathfrak f 0} +  I_{\mathfrak f 1} N_{\mathfrak f 0} \\
\psi_{\mathfrak f b} &= I_{\mathfrak f 0} \psi_{\mathfrak f 0} + I_{\mathfrak f 1} N_{\mathfrak f a} \\
\psi_{\mathfrak f c} &= I_{\mathfrak f 0} \psi_{\mathfrak f a} + I_{\mathfrak f 1} (2 N_{\mathfrak f b} - N_{\mathfrak f 0}) \\
\psi_{\mathfrak f} (t+\Delta t) &= I_{\mathfrak f 2} \psi_{\mathfrak f 0} + I_{\mathfrak f 3} N_{\mathfrak f 0} + I_{\mathfrak f 4} (N_{\mathfrak f a} + N_{\mathfrak f b}) + I_{\mathfrak f 5} N_{\mathfrak f c}
\end{aligned}
$$

where

$$
\begin{aligned}
I_{\mathfrak f 0} &= e^{\omega_{\mathfrak f} \Delta t/2} \\
I_{\mathfrak f 1} &= \frac{1}{\omega_{\mathfrak f}}
( e^{ \omega_{\mathfrak f} \Delta t/2} - 1) \\
I_{\mathfrak f 2} &= e^{\omega_{\mathfrak f} \Delta t} \\
I_{\mathfrak f 3} &= \frac{1}{ \omega_{\mathfrak f}^3\Delta t^2}
\left ( -4 -  \omega_{\mathfrak f} \Delta t  + e^{\omega_{\mathfrak f} \Delta t}(4-3\omega_{\mathfrak f} \Delta t + \omega_{\mathfrak f}^2 \Delta t^2 )  \right ) \\
I_{\mathfrak f 4} &= \frac{2}{ \omega_{\mathfrak f}^3\Delta t^2}
\left ( 2 + \omega_{\mathfrak f} \Delta t + e^{\omega_{\mathfrak f} \Delta t}(-2 + \omega_{\mathfrak f} \Delta t) \right ) \\
I_{\mathfrak f 5} &= \frac{1}{ \omega_{\mathfrak f}^3\Delta t^2}
\left ( -4 - 3 \omega_{\mathfrak f} \Delta t -  \omega_{\mathfrak f}^2 \Delta t^2 + e^{\omega_{\mathfrak f} \Delta t}(4-\omega_{\mathfrak f} \Delta t) \right )
\end{aligned}
$$

---

**Algorithm:** The ETD4RK scheme

In the small $\omega_{\mathfrak f}$ limit, we have

$$
I_{\mathfrak f 0} \approx 1
$$

$$
I_{\mathfrak f 1} \approx \frac{1}{2} \Delta t
$$

$$
I_{\mathfrak f 2} \approx 1
$$

$$
I_{\mathfrak f 3} \approx
\frac{1}{ \omega_{\mathfrak f}^3\Delta t^2} \times
\left ( -4 - \omega_{\mathfrak f} \Delta t + (1 + \omega_{\mathfrak f} \Delta t + \frac{1}{2} (\omega_{\mathfrak f} \Delta t)^2 + \frac{1}{6} (\omega_{\mathfrak f} \Delta t)^3 )(4-3\omega_{\mathfrak f} \Delta t + \omega_{\mathfrak f}^2 \Delta t^2 )
\right ) $$

$$
= \frac{1}{ \omega_{\mathfrak f}^3\Delta t^2}
\left ( \frac{4}{6} (\omega_{\mathfrak f} \Delta t)^3 - \frac{3}{2} (\omega_{\mathfrak f} \Delta t)^3 + (\omega_{\mathfrak f} \Delta t)^3
\right )
= \frac{1}{6} \Delta t
$$

$$
I_{\mathfrak f 4} \approx \frac{2}{ \omega_{\mathfrak f}^3\Delta t^2}
\left (2 + \omega_{\mathfrak f} \Delta t +(1 + \omega_{\mathfrak f} \Delta t + \frac{1}{2} (\omega_{\mathfrak f} \Delta t)^2 + \frac{1}{6} (\omega_{\mathfrak f} \Delta t)^3 )(-2 + \omega_{\mathfrak f} \Delta t)
\right )
$$

$$
= \frac{2}{ \omega_{\mathfrak f}^3\Delta t^2}
\left ( \frac{1}{2} (\omega_{\mathfrak f} \Delta t)^3-\frac{2}{6}(\omega_{\mathfrak f} \Delta t)^3
\right ) = \frac{1}{3} \Delta t
$$

$$
I_{\mathfrak f 5} =
\frac{1}{ \omega_{\mathfrak f}^3\Delta t^2} \times
\left (
-4 - 3 \omega_{\mathfrak f} \Delta t -  \omega_{\mathfrak f}^2 \Delta t^2 + (1 + \omega_{\mathfrak f} \Delta t + \frac{1}{2} (\omega_{\mathfrak f} \Delta t)^2 + \frac{1}{6} (\omega_{\mathfrak f} \Delta t)^3 )(4-\omega_{\mathfrak f} \Delta t)
\right )  
$$

$$
=\frac{1}{ \omega_{\mathfrak f}^3\Delta t^2}
\left ( \frac{4}{6} (\omega_{\mathfrak f} \Delta t)^3 - \frac{1}{2} (\omega_{\mathfrak f} \Delta t)^3
\right ) = \frac{1}{6} \Delta t
$$

Similar as for the EDT2RK case $I_{\mathfrak f 1}$, $I_{\mathfrak f 3}$, $I_{\mathfrak f 4}$, and $I_{\mathfrak f 5}$ contains a division by $0$ when $\omega_{\mathfrak f} = 0$.  
We therfore replace these coeficients with their limits when $|\omega_{\mathfrak f}|$ is smaller than a tolerance.
This has been important in order to make the the code stable for some of the systems.
In the same way as the EDT2RK scheme this is implemented as the function
`self.evolve_ETD4RK_loop(self, integrating_factors_f, non_linear_evolution_function, field, field_f)`
This function is called by the evolvers discussed in the model chapter if the method is defined as ```method = "ETD4RK"```, the integrating factors are found with
`self.calc_evolution_integrating_factors_ETD4RK(self, omega_f, tol=10**(-4))`.

### The fully non-linear limit
It is both interesting and enlightening to see the fully non-linear limit of these equations, i.e., the limit in which $\omega_{\mathfrak f} =0$, $N_{\mathfrak f} = \partial_t \psi \equiv \dot{\psi}_{\mathfrak f}$ and the small $\omega_{\mathfrak f}$ approximations are exact.
For the ETD2RK scheme, we get

$$
\psi_{\mathfrak f a} = \psi_{\mathfrak f 0} +  \dot{\psi}_{0_f} \Delta t
$$

$$
\psi(t+\Delta t)_f = \psi_{\mathfrak f 0} + \dot{\psi}_{\mathfrak f 0} \frac{\Delta t}{2} + \dot{\psi}_{\mathfrak f a} \frac{\Delta t}{2},
$$

which is a two-stage Runge-Kutta method called Heun's method.

The ETD4RK scheme becomes

$$
\psi_{\mathfrak f a} = \psi_{\mathfrak f 0} +  \dot{\psi}_{\mathfrak f 0} \frac{\Delta t}{2}
$$

$$
\psi_{\mathfrak f b} =  \psi_{\mathfrak f 0} +  \dot{\psi}_{\mathfrak f a} \frac{\Delta t}{2}
$$

$$
\psi_{\mathfrak f c} = \psi_{\mathfrak f a} + ( 2 \dot{\psi}_{\mathfrak f b} - \dot{\psi}_{\mathfrak f 0}) \frac{\Delta t}{2}
$$

$$
\psi_{\mathfrak f} (t+\Delta t) = \psi_{\mathfrak f 0} + \frac{1}{6} ( \dot{\psi}_{\mathfrak f 0} + 2 \dot{\psi}_{\mathfrak f a} + 2 \dot{\psi}_{\mathfrak f b} + \dot{\psi}_{\mathfrak f c} ) \Delta t.
$$

Note that this is not the typical Runge-Kutta 4 method, due to the differences in calculating $\psi_{\mathfrak f c}$.
The reason is that a straight-forward generalization of the Runge-Kutta 4 method will not produce a fourth-order method in the general case [^coxExponentialTimeDifferencing2002].

### The fully linear limit

If $N=0$, the evolution equation changes to

$$
\psi_{\mathfrak f}(t+\Delta t) = e^{\omega_{\mathfrak f} \Delta t} \psi_{\mathfrak f}.
$$

An example of this is the schrödinger equation, for which $\omega_{\mathfrak f} = -\frac{1}{2}  i \mathbf k^2$, so we get

$$
\psi_{\mathfrak f}(t+\Delta t) = e^{- i \frac{1}{2} \mathbf k^2 \Delta t} \psi_{\mathfrak f}.
$$

This is an exact eqution, of course, so you may evolve this free particle solution to any time.

## Testing

In order to test the numerical methods, study the simplest model of a field equation with a (non-linear) forcing term, namely the heat equation

$$
\partial_t T = \nabla^2 T + f(\mathbf r),
$$

where $T$ is the temperature in celcius, and $f(\mathbf r)$ is a forcing term, which we model as

$$
f(\mathbf r) = A (T_0-T) \exp\left (-\frac{(\mathbf r-\mathbf r_0)^2}{2\sigma^2}\right ),
$$

which represents a heating element with temperature $T_0$ placed at $\mathbf r_0$.

As a benchmark, we use the `solve_ivp` of the scipy library `sp.integrate` to solve the equation using a finite difference method.
The solutions match to a satisfactory degree, but a more thorough investigation into how the accuracy of the framework and integration methods scale with spatial and temporal resolution will be performed in the future.
Tests are included in `test_base_system.py`, but for visual examination, here are animations of the initial condition $T=0$ in all three dimensions

![Testing of the evolution code in 1 dimension](images/base_system_evolution_test_1D.gif#only-light)
![Testing of the evolution code in 1 dimension](images/base_system_evolution_test_1D-colorinverted.gif#only-dark)

![Testing of the evolution code in 2 dimensions](images/base_system_evolution_test_2D.gif#only-light)
![Testing of the evolution code in 1 dimension](images/base_system_evolution_test_2D-colorinverted.gif#only-dark)

![Testing of the evolution code in 3 dimensions](images/base_system_evolution_test_3D.gif#only-light)
![Testing of the evolution code in 1 dimension](images/base_system_evolution_test_3D-colorinverted.gif#only-dark)

## Algorithms for tracking defects

To be written

## Calculating the velocity

The equations for the velocity are taken from Ref.[^skogvollUnifiedFieldTheory2023], simplified using Mathematica and then
substituted for python code using chatGPT.


## Tutorial: Basic framework

This tutorial is meant to introduce the very basic functionality of the comfit Python package.

The class underlying everything is the class called `BaseSystem`. Every other model (classes located in the `./models/` folder) are subclasses of BaseSystem.

Let's start by importing the package

```python
import comfit as cf
```

Now we can define a system with a given dimension and resolution

```python
sys1 = cf.BaseSystem(1,xRes=11)
```

`sys1` now contains a number of useful parameters and functions, lets print some of them.

```python
print(sys1.x)
```

```python
print(sys1.dx, sys1.dt)
```

Even though the system is one-dimensional in this case, it still contains some dummy variables related to y and z.

```python
print(sys1.y,sys1.dy,sys1.z,sys1.dz)
```

Of particular interest, perhaps, is the length scale parameter 

```python
print(sys1.a0)
```

Now, we can try to do some basic calculating and plotting

```python
import numpy as np
y = np.sin(sys1.x)
sys1.plot_field(y)
```

By changing the length scale to $a_0=2\pi$, we get something that clearly shows the periodicity.

```pyhthon
sys1.a0 = 2* np.pi
sys1.plot_field(y)
```

The function looks quite jagged, which we can fix by increasing the resolution and decreasing the interspacing dx

```python
sys2 = cf.BaseSystem(1,xRes=101,dx=0.1)
y2 = np.sin(sys2.x)
sys2.plot_field(y2)
```

Lets have a look at a 2D system and plot a 2D function

```py
sys2 = cf.BaseSystem(2,xRes=101,dx=0.1,yRes=101,dy=0.1)
field = sys2.x-1 + 0.7*sys2.y-5
sys2.plot_field(field)
```

Now, the very attentive python reader will have noticed something strange about the previous cell. Namely how is it possible that by adding two 1D arrays 
 and 
 do we produce a 2D array. The key lies in the way that these are stored. Let's look at them.

```python
print("x:\n",sys2.x)
```
```python
print("y:\n",sys2.y)
```
As you see, the x-array is stored as a 
 array, while the y-array is stored as 
. When doing calculations with arrays in this way, it is not necessary to create a meshgrid to produce a 2D matrix.

### Exercise 1
Plot the function

$$
f(x,y) = \sin(x+y)
$$

for $x \in [-1,1]$, $y \in [-1,1]$.

??? note "Solution"
    ```python
    import comfit as cf
    import numpy as np
    import matplotlib.pyplot as plt

    bs = cf.BaseSystem(2, xlim=[-3,3], ylim=[-3,3])
    f = np.sin(bs.x+bs.y)
    bs.plot_field(f)
    plt.show()
    ```

### Exercise 2
Plot the function

$$
f(x,y,z) = \exp(-(x^2 + y^2 + z^2))
$$

for $x \in [-3,3]$, $y \in [-3,3]$ and $z \in [-3,3]$. 

??? note "Solution"
    ```python
    import comfit as cf
    import numpy as np
    import matplotlib.pyplot as plt

    bs = cf.BaseSystem(3, xlim=[-3,3], ylim=[-3,3], zlim=[-3,3])
    f = np.exp(-(bs.x**2+bs.y**2+bs.z**2))
    bs.plot_field(f,number_of_layers=3)
    plt.show()
    ```

Now, plot the funciton in the plane through $(0,0,0)$ given by the normal vector $[1,1,1]$ 

??? note "Solution"
    ```python
    import comfit as cf
    import numpy as np
    import matplotlib.pyplot as plt

    bs = cf.BaseSystem(3, xlim=[-3,3], ylim=[-3,3], zlim=[-3,3])
    f = np.exp(-(bs.x**2+bs.y**2+bs.z**2))
    bs.plot_field_in_plane(f, position=[0,0,0], normal_vector=[1,1,1],plotting_lib='matplotlib')
    plt.show()
    ```

[^coxExponentialTimeDifferencing2002]: Cox, S. M., & Matthews, P. C. (2002). Exponential Time Differencing for Stiff Systems. Journal of Computational Physics, 176(2), 430–455. [https://doi.org/10.1006/jcph.2002.6995](https://doi.org/10.1006/jcph.2002.6995)

[^skogvollUnifiedFieldTheory2023]: Skogvoll, V., Rønning, J., Salvalaglio, M., & Angheluta, L. (2023). A unified field theory of topological defects and non-linear local excitations. Npj Computational Materials, 9(1), Article 1. [https://doi.org/10.1038/s41524-023-01077-6](https://doi.org/10.1038/s41524-023-01077-6)

[^skogvollPhaseFieldCrystal2022]: Skogvoll, V., Angheluta, L., Skaugen, A., Salvalaglio, M., & Viñals, J. (2022). A phase field crystal theory of the kinematics of dislocation lines. Journal of the Mechanics and Physics of Solids, 166, 104932. https://doi.org/10.1016/j.jmps.2022.104932
