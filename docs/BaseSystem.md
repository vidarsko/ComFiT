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

where $r^2 = (x-\xmid)^2 + (y-\ymid)^2$, which is a function that is zero in the center region and goes to $1$ at infinity over a width of $w$. The filtered field $\eta$ is then obtained by making a smooth function that goes from the previous angle field according to

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

where $\omega$ is a linear differential operator and $N$ is a non-linear operator (function of $\psi$). Table (referenced as `tab:time_evo_operator_examples_non_conserved`) shows some examples from the models that we will discuss in the following chapters.

| Model | $\omega$ | $\omega_f(\mathbf{k})$ | $N$ |
| --- | --- | --- | --- |
| Quantum Mechanics | $\frac{1}{2}i \nabla^2 $ | $-\frac{1}{2} i \mathbf{k}^2$ | $- i V$ |
| BEC | $(i+\gamma) (1+\frac{1}{2}\nabla^2)$ | $(i+\gamma) (1-\frac{1}{2}k^2)$ | $- (i + \gamma) (V_{ext} + |\psi|^2)\psi$ |
| Active Nematic | $\frac{K}{\gamma} \nabla^2 +\frac{AB}{\gamma}$ | $-\frac{K}{\gamma} k^2 +\frac{AB}{\gamma}$ | $- \mathbf{u}\cdot \nabla Q + Q \Omega -\Omega Q - \frac{2A}{\gamma}Q^2_{kk}Q$ |

*Table: Examples of time evolution operators, non-conserved.*

In the following, we will explain the method of evolution of exponential time differencing as introduced in Ref. [coxExponentialTimeDifferencing2002] for stiff systems. This will result in two integration schemes, the exponential time differencing second order Runge Kutta 2 (ETD2RK) scheme and the forth order ETD4RK scheme. As in Ref. [coxExponentialTimeDifferencing2002], we will show an intuitive way to obtain the former and only recite the expressions for the latter.

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
\psi_f (t+\Delta t) = \psi_f (t) e^{ \omega_f \Delta t} + e^{\omega_f \Delta t} \frac{1}{- \omega_f} [e^{- \omega_f \tau}]_0^{\Delta t} N_{f0} + e^{ \omega_f \Delta t} \frac{1}{\Delta t} [\frac{\tau e^{-\omega_f \tau}}{-\omega_f} - \frac{e^{-\omega_f \tau}}{\omega_f^2}]_0^{\Delta t} \Delta N_f
$$

To find $\psi_f (t+\Delta t)$, we would need to know the value $N_f (t+\Delta t)$ before finding the state at $\psi(t+\Delta t)$. To do this, we first find a predicted state $\psi_a$ by assuming $\Delta N_f=0$ and calculating $\psi(t)$ according to the equation above. This lets us calculate an approximate $\Delta N_f = N_{fa} - N_{f0}$ and we use this in order to evolve $\psi$. This is the ETD2RK scheme.

$$
\psi_{fa} = \psi_{f0} + \psi_{f1} \psi_{f0}
$$

$$
\psi_f (t+\Delta t) = \psi_{fa} + \psi_{f2} (N_{fa} - N_{f0})
$$

where

$$
\psi_{f0} = e^{\omega_f \Delta t}
$$

$$
\psi_{f1} = \frac{1}{\omega_f} (e^{ \omega_f \Delta t} - 1)
$$

$$
\psi_{f2} = \frac{1}{\Delta t \omega_f^2} (e^{ \omega_f \Delta t} -1  -\omega_f \Delta t)
$$

For numerical purposes, it is useful to calculate the small $\omega_f$ limit. We expand the exponential in its Taylor series and keep the leading order term to get:

$$
\psi_{f0} \approx 1
$$

$$
\psi_{f1} \approx \Delta t
$$

$$
\psi_{f2} \approx \frac{1}{2} \Delta t
$$

In $\psi_{f1}$, and $\psi_{f2}$ there is a division by $0$ when $\omega_f = 0$. To avoid numerical issues related to this we use the above limits when $|\omega_f|$ is smaller than a tolerance. We don't use the limit for $\psi_{f0}$ since it doesn't contain a division by $0$. The function `evolve_ETD2RK_loop` defined in the base system class performs an ETD2RK step. This function is called by the evolvers discussed in the model chapter if the method is defined as `method = "ETD2RK"`. This is the default solver if `method` is not set. The integrating factors for a given $\omega_f(\mathbf{k})$ can be found with the function `calc_evolution_integrating_factors_ETD2RK` where the variable `tol` gives when the factors should be replaced by their leading order Taylor expansion.

(Note: The content is lengthy and includes complex mathematical expressions. It might not render correctly on all Markdown platforms, especially those that do not support LaTeX math formatting. The code snippets and tables are formatted as standard Markdown.)


# Algorithms for tracking defects {#Chapter:defect_tracking}

## Tracking the indices?

There is a question to be settled as to how we should track the zeros,
because

## Methods for tracking

There are multiple possible approaches to tracking defects. This is
because we have access to both a coarse defect field, a singular defect
field and the zeros of the wave function.

### Method 1: Using the zeros of the wave function and the sign of the coarse defect density

In this case, we consider the zeros of the order parameter as the
indicator of the location of the defects and use the value of the defect
density to calculate the nature of the defect.

Advantages: An advantage of this method is that it captures the nature
of defects very close up to the point of annihilation. Also, it provides
a very natural way to find the position of the defect without being
fixed to the predefined coordinates.

### Method 2: Integrating the singular defect density field

Advantage: This method will allow for an intuitive emptying of the
field, adding defects one by one untill they annihilate.

Disadvantage: Due to the calculation of the delta function, the defect
density becomes a very sharply peaked function that might land inbetween
position indices. Not optimal

### Method 3: Integrating the defect density field

Advantage: As long as the defects are far apart, this method will give a
nice and concrete way of calculating them.

Disadvantage: In order to get the full charge of the defect, you will
need to integrate over quite a large region, which is not optimal. Also,
as the defects get close, they will start to overlap, giving defects
with partial charges.

From me messing around with the code, it actually seems as this is the
most fruitful approach. I think we need to quantify the value at which
the defects are no longer isolated, i.e., how large of a radius we need
to include in order to capture the charge of the defect.

## real method

For all the systems, we calcuate a quanity which when integrated gives
the charge.

After integrating, if the value exceeds the threshold, we remove the
area. We also limit the search region to outside a ball of 2 rad, to
avoid halting the search algorithm.

in the case of 3d, the charge will scale with a0 squared

## Calculating the velocity

The equations for the velocity are taken from Ref.
[@skogvollUnifiedFieldTheory2023], simplified using Mathematica and then
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
instances, line2D objects or text. Ref. [@yeManyWaysCall2020] gives a
pedagogical overview, and Fig.
[9.1](#fig:MPLFigureAndAxes){reference-type="ref"
reference="fig:MPLFigureAndAxes"} shows the idea.

![The different levels of plotting in matplotlib. Figure reprinted from
Ref. [@yeManyWaysCall2020] with permission (not yet confirmed).
](Figures/FigureAndAxes.png){#fig:MPLFigureAndAxes width="\\textwidth"}

Therefore, the implemented plot functions will all take the ax as an
optional input and give that as an optional output if not provided.

Each class has a function called `plot`, which is our best idea of how
to plot the current configuration of the field.

## Predefined figure sizes

Obviously, a lot of these plots are meant for publication. Therefore,
there is a tool package that contains predefined and fitting sizes for
the figures in question. \[TO BE ADDED\]

## The plot in plane function

The following algorithm is implemented for the plot in plane function.
The goal is to plot the field on a plane $P$ perpendicular to a normal
vector $\boldsymbol{n}$ through a point $P_0$.

The fundamental problem of the plot in plane function is that the field
is only defined on a grid, a
$(\textrm{\lstinline|xRes|},\textrm{\lstinline|yRes|},\textrm{\lstinline|zRes|})$
matrix `field`.

Given the normal vector to the plane $P$ on which we want to plot the
function, we create a grid of x- and y-coordinates $X_0$ and $Y_0$ and
generate a vector $\boldsymbol{R}$, which is the x- y- and z-position of
the points on the grid. In order to make sure that we cover the whole
domain, we first generate points of $\boldsymbol{R}$ that are also
outside of the domain and then remove those excess points outside of the
domain.

Thus, all the coordinates in $\boldsymbol{R}$ are inside the domain
covered by the matrix. Given such a point $\boldsymbol{R}$, we find the
float index `Ri` corresponding to the point, e.g., if $R_x=1.7$ and we
have `x[2] = 1` and `x[3] = 2`, then $\textrm{\lstinline|Ri[0]|} = 2.7$.
From this vector, we determine the indices of the corners `RiLH` of the
voxel containing $\boldsymbol{R}$ and the (positive) index distances
`dRi` from these corners to `Ri`, see Fig.
[\[fig:PlotInPlaneFigure\]](#fig:PlotInPlaneFigure){reference-type="ref"
reference="fig:PlotInPlaneFigure"}.

::: overpic
Figures/PlotInPlaneIllustration/PlotInPlaneIllustration.png
(0,3)(`RiLH[0][0]`,`RiLH[1][0]`,`RiLH[2][0]`)
(65,3)(`RiLH[0][1]`,`RiLH[1][0]`,`RiLH[2][0]`) (23,13)`dRi[0][0]`
(49,13)`dRi[0][1]` (66,17)`dRi[1][0]` (81,24)`dRi[1][1]`
(-3,34)`dRi[2][0]` (-3,58)`dRi[2][0]` (38,55)(`Ri[0]`,`Ri[2]`,`Ri[3]`)
:::

We then calculate the value `field_on_plane` to be plotted at
$\boldsymbol{R}$ by a weighted sum $$\begin{gathered}
\textrm{\lstinline{field_on_plane}} \\
=
\sum_{\textrm{\lstinline|ix|,\lstinline|iy|,\lstinline|iz|}=0}^1
(1 - \textrm{\lstinline|dRi[0][ix]|})
(1 - \textrm{\lstinline|dRi[1][iy]|})
(1 - \textrm{\lstinline|dRi[2][iz]|})
\textrm{\lstinline|field[ix,iy,iz]|},
\end{gathered}$$ which corresponds to a sum of the field at the voxel
corners weighted by the index distance to `Ri`.

## Animation

Easy animation

# How to create your own model

Creating your own class with this framework is easy. Here is a
step-by-step instruction on how to do it.

::: tcolorbox
Example system: The Landau theory of the Bragg-Williams theory
[@chaikinPrinciplesCondensedMatter1995]
$$\mathcal F[\phi] = \int d\boldsymbol{r} \frac{1}{2} \texttt r \phi^2 - \texttt w \phi^3 + \texttt u \phi^4 + \frac{1}{2} \texttt c (\nabla \phi)^2.$$
Seeking to minimize this functional in equilibrium, we have a simple
equation of motion
$$\partial_t \phi = - \nabla^2 \frac{\delta \mathcal F}{\delta \phi}.$$
So, before you start, you need to create a class structure for your
system which inherits the basesystem class and sets these parameters. It
could look something like this.

    import comfit as cf

    class Landau(BaseSystem):
    """
    Class that ...
    """
    def __init__():
        something
:::

## Step 1: Rephrase your problem

Write the problem on the form $\partial_t \psi = \omega \psi +N(\psi)$,
and $N$ might be dependent on other fields as well.

$\psi$ may well be a multi-component quantity.

::: tcolorbox
In the previous example, we had the equation of motion. Expanding that
we get $$\partial_t \phi = ...$$ We see that we have a linear part, and
a non-linear part $$\omega \phi =$$ $$N =$$
:::

## Step 2: 

Put the expression for ${{\omega}_{\scriptscriptstyle \mathbbm f}}$ into
the function that calculates the integrating factors for EDT2RK or
EDT4RK.

::: tcolorbox
It could be

    def calc_omega_f(self):
        Something
:::

## Step 3: Create the non-linear function

Make a function that calculates the non-linear function $N(\psi)$ in
Fourier space.

::: tcolorbox
In this case, it could look like

    def calc_non_linear_part(self,field):
        Something
:::

## Step 4: Make the evolver

Make an evolver by inserting the integrating factors and the non-linear
evolution function into the evolver loop corresponding to the
integrating factors you found.

    def evolver(self, number_of_steps, method='ETD2RK'):

        omega_f = ...

        integrating_factors_f, solver = self.calc_integrating_factors_f_and_solver(omega_f, method)

        for n in tqdm(range(number_of_steps), desc='Evolving the model'):
            self.psi, self.psi_f = solver(integrating_factors_f,
                                          self.calc_nonlinear_evolution_function_f,
                                          self.psi, self.psi_f)
            self.psi = np.real(self.psi) #If your field is real

::: tcolorbox
So in this case, it would be

    def evolve_landau(self,number_of_steps):
        solver = ...
:::

## Step 5: Configure some initial condition

::: tcolorbox
So in this case, it would be

    def conf_initial_condition(self,...):
        self.phi = np.zeros()
:::
