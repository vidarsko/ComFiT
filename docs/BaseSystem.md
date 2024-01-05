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
| $\vec r_{\textrm{mid}}$ | `rmid`      | |
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
\int d^d r' \mathcal K(\vec r-\vec r') 
\tilde \rho(\vec r') ,
$$

where $\mathcal K(\vec r'-\vec r)$ is a Gaussian kernel given by 
From a numerical point of view, this is done in Fourier space since, by the convolution theorem, 

$$
\mathcal F \rho = \mathcal F {\mathcal K} \mathcal F {\tilde \rho}.
$$

Thus, we need the Fourier transform of $\mathcal K$, which is 

$$
\mathcal F {\mathcal K} = \int d^d r e^{-i \vec k \cdot \vec r} \frac{1}{(2\pi w^2)^{d/2}} \exp\left (-\frac{\vec r^2}{2w^2} \right )
= \frac{1}{(2\pi w^2)^{d/2}} \prod_{n=1}^d \int dr_n e^{-\frac{1}{2 w^2} r_n^2 -i k_n r_n} 
= \frac{1}{(2\pi w^2)^{d/2}} \prod_{n=1}^d \int dr_n e^{-\frac{1}{2 w^2} (r_n^2 + 2 i w^2 k_n r_n)}  
= \frac{1}{(2\pi w^2)^{d/2}} \prod_{n=1}^d e^{-\frac{1}{2} w^2 k_n^2} \int dr_n e^{-\frac{1}{2 w^2} (r_n + i w^2 k_n)^2}  
= e^{-\frac{1}{2} w^2 \vec k^2}.
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

In two dimensions, the angle field takes values $\theta \in [-\pi,\pi \rangle$ and a vortex is a point $\vec r_0$.
The angle field produced by the vortex has a circulation which is a multiple integer of $2\pi$, i.e.

$$
\oint d\theta = 2\pi s_n, 
$$

where $s_n$ is the integer charge of the vortex. 
A possible angle field for a vortex positioned at $(x_0,y_0)$ is given by 

$$
\theta_n = s_n \textrm{atan2}(y-y_0,x-x_0)
$$
