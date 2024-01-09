# Folders and file types

* Documentation is written in markdown.
* Tutorials are written in markdown.

# File naming conventions

* Documentation files are named using CamelCase. 
* Image files in the docs folder are named using `snake_case`

# Mathematical conventions



The imaginary unit will be denoted $\mathfrak i$, `\mathfrak i` to avoid confusion with the index $i$. 
Index symmetrization $()$ and anti-symmetrization $[]$ will be used throughout.
They are defined for an tensor $A$ with two indices by 

$$
A_{(ij)} = \frac{1}{2} (A_{ij} + A_{ji}),
$$

$$
A_{[ij]} = \frac{1}{2} (A_{ij} - A_{ji})
$$


The Fourier transform of a field $\psi$ will be denoted $\psi_{\mathfrak f}$, `\psi_{\mathfrak f}`, and is defined as 

$$
\psi_{\mathfrak f} (\mathbf k) = \int d^d r e^{-\mathfrak i \mathbf k \cdot \mathbf r} f(\mathbf r),
$$

and the inverse is given by 

$$
\psi(\mathbf r) = \frac{1}{(2\pi)^d} \int d^d k e^{\mathfrak i\mathbf k\cdot \mathbf r} \psi_{\mathfrak f}(\mathbf k),
$$

where $d$ is the spatial dimension.

Vectors $\mathbf a, \mathbf b, \mathbf c, \boldsymbol \Omega$ are denoted using boldfont (`\mathbf`, `\boldsymbol`) , while rank 2 tensors vary more. 
Typical choices however are non-bold greek letters ($\sigma$) lower case Fraktur letters ($\mathfrak h$, `\mathfrak h`) or capital letters ($Q$).

The dot product ($\cdot$) is a contraction over the last index

$$
(\nabla \cdot \sigma)_i = \partial_j \sigma_{ij},
$$

while the double dot product $\dot \cdot$ is a contraction over the last two indices

$$
(\mathcal C \cdot \cdot \mathfrak e)_{ij} = \mathcal C_{ijkl} \mathfrak e_{kl}
$$

$$
\psi(\mathbf r) = \frac{1}{(2\pi)^d} \int d^d k e^{\mathfrak i\mathbf k\cdot \mathbf r} \psi_{\mathfrak f}(\mathbf k),
$$


# Programming notation conventions

*  [PEP8](https://peps.python.org/pep-0008/) for python programming, and [PEP257](https://peps.python.org/pep-0257/)/[Google Python Style Guide](https://github.com/google/styleguide/blob/gh-pages/pyguide.md) for doc strings.
