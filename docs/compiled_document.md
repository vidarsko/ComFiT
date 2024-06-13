# Compiled Markdown Document

## Unnamed Section

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


