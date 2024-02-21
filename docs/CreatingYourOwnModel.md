# Creating your own model

Using ComFiT to create your own model is easy.

## Level 1: Get it up and running
Creating your own class with this framework is easy. Here is a
step-by-step instruction on how to do it.

Consider as an instructive example, the Landau theory of the Bragg-Williams theory [^chaikinPrinciplesCondensedMatter1995]

<div style="border:2px solid black; padding:10px; margin:10px;">

$$
\mathcal F[\phi] = \int d\boldsymbol{r} \frac{1}{2} \texttt r \phi^2 - \texttt w \phi^3 + \texttt u \phi^4 + \frac{1}{2} \texttt c (\nabla \phi)^2.
$$

Seeking to minimize this functional in equilibrium, we have a simple
equation of motion
$$\partial_t \phi = - \nabla^2 \frac{\delta \mathcal F}{\delta \phi}.$$
So, before you start, you need to create a class structure for your
system which inherits the basesystem class and sets these parameters. It
could look something like this.
```python
import comfit as cf

class Landau(BaseSystem):
"""
Class that ...
"""
def __init__():
    something
```
</div>


## Step 1: Rephrase your problem

Write the problem on the form $\partial_t \psi = \omega \psi +N(\psi)$,
and $N$ might be dependent on other fields as well.

$\psi$ may well be a multi-component quantity.

<div style="border:1px solid black; padding:10px;">


In the previous example, we had the equation of motion. Expanding that
we get $$\partial_t \phi = ...$$ We see that we have a linear part, and
a non-linear part $$\omega \phi =$$ $$N =$$
</div>

## Step 2: 

Put the expression for ${{\omega}_{f}}$ into
the function that calculates the integrating factors for EDT2RK or
EDT4RK.

<div style="border:1px solid black; padding:10px;">
It could be

```python
def calc_omega_f(self):
        Something
```
    
</div>


## Step 3: Create the non-linear function

Make a function that calculates the non-linear function $N(\psi)$ in
Fourier space.

<div style="border:1px solid black; padding:10px;">
In this case, it could look like

```python
def calc_nonlinear_part(self,field,t):
    Something
```
</div>

## Step 4: Make the evolver

Make an evolver by inserting the integrating factors and the non-linear evolution function into the evolver loop corresponding to the integrating factors you found.

    def evolver(self, number_of_steps, method='ETD2RK'):

        omega_f = ...

        integrating_factors_f, solver = self.calc_integrating_factors_{fa}nd_solver(omega_f, method)

        for n in tqdm(range(number_of_steps), desc='Evolving the model'):
            self.psi, self.psi_f = solver(integrating_factors_f,
                                          self.calc_nonlinear_evolution_function_f,
                                          self.psi, self.psi_f)
            self.psi = np.real(self.psi) #If your field is real

<div style="border:2px solid black; padding:10px; margin:10px;">
So in this case, it would be

    def evolve_landau(self,number_of_steps):
        solver = ...
</div>


## Step 5: Configure some initial condition

<div style="border:1px solid black; padding:10px;">

```python
def conf_initial_condition(self,...):
        self.phi = np.zeros()
```
Your content goes here.
</div>

## Level 2: Merge it with the source


If you want to contribute to the source code with your model (which is highly appreciated) to the models folder, you need to do the following while you are developing your code. Follow the examples provided by the `core_models`.

* Document while writing the code, following the [notation conventions](Conventions.md).
* Write documentation for the code, explaining the mathematics of the model. Remember to follow the [mathematical notation convention](Conventions.md#mathematical-notation-convention).
* Make tests for the model.
* [optional] Create tutorials showcasing how to use you code.


[^chaikinPrinciplesCondensedMatter1995]: 