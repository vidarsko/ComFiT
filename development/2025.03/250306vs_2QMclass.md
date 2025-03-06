# A 2-QM particle class

This is a document for me to see how it would be if I were to simulate a 2-QM particle class.

In this case, the wavefunction would be a funciton of two variables:

$$
\psi(x_1, x_2)
$$

The Hamiltonian would be a function of the two variables:

$$
H(x_1, x_2) = -\frac{\hbar^2}{2m_1}\nabla_1^2 -\frac{\hbar^2}{2m_2}\nabla_2^2 + V(x_1, x_2)
$$

if they are non-interacting, then the potential would be a sum of the two potentials:

$$
V(x_1, x_2) = V_1(x_1) + V_2(x_2)
$$

So this is what I need to implement.

$$
\partial_t \psi(x_1, x_2) = \mathfrak i 
\left (\frac{1}{2}\nabla_1^2 + \frac{1}{2} \nabla_2^2 
 - V(x_1,x_2)\right )\psi(x_1, x_2)
$$