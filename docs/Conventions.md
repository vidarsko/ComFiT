# Conventions

In this section, we will describe the conventions used in the documentation and code.

## Folders and file types

* Documentation is written in markdown.
* Tutorials are written in markdown.

## File naming conventions

* Documentation files are named using CamelCase.
* Image files in the docs folder are named using `snake_case`

## Mathematical conventions

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

```math
(\nabla \cdot \sigma)_i = \partial_j {\sigma}_{ij}
```

while the double dot product $\dot \cdot$ is a contraction over the last two indices

$$
(\mathcal C \dot \cdot \mathfrak e)_ {ij} = \mathcal C_ {ijkl} \mathfrak e_ {kl}
$$

## Programming notation conventions

* [PEP8](https://peps.python.org/pep-0008/) for python programming
* [NumPy docstring format](https://numpydoc.readthedocs.io/en/latest/format.html) for documentation strings
* Import ordering should follow this pattern:
  1. Standard library imports
  2. Third-party library imports
  3. Local application imports
  4. Each group should be separated by a blank line

Stand-alone functions are documented as follows:

```python
def function_name(
        arg1: arg1_type, 
        arg2: arg2_type, 
        arg3: Optional[arg3_type] = None,
        **kwargs: Any
        ) -> return_type:
    """Short description of the function (in the imperative mood).

    Optional longer description of the function.

    Parameters
    ----------
    arg1 : arg1_type
        Description of arg1. No need to state default values.
    arg2 : arg2_type
        Description of arg2.
    arg3 : arg3_type, optional
        Description of arg3. Defaults to None.
    \*\*kwargs : Any
        Description of additional keyword arguments.
    

    Returns
    -------
    return_type
        Description of the return value.

    Raises
    ------
    Exception
        Description of the exception.

    Examples
    --------
    >>> function_name(1, 2)
    3       
    """
    pass
```

If the function is a method of a `comfit` model, it should be documented as follows:

```python
from typing import TYPE_CHECKING # Import necessary typing packages

if TYPE_CHECKING:
    from comfit.core.base_system import BaseSystem

# General packages
# Import necessary packages from the standard library, e.g.
# import numpy as np

# Comfit packages
# Import necessary packages from comfit from the subpackages, e.g. 
# from comfit.core import BaseSystem

def function_name(
        self: 'BaseSystem',
        arg1: arg1_type, 
        arg2: arg2_type, 
        arg3: Optional[arg3_type] = None
        ) -> return_type:
    """Short description of the function (in the imperative mood).

    Optional longer description of the function.

    Parameters
    ----------
    arg1 : arg1_type
        Description of arg1. No need to state default values.
    arg2 : arg2_type
        Description of arg2.
    arg3 : arg3_type, optional
        Description of arg3. Defaults to None.
    \*\*kwargs : Any (use backslash to escape the asterisk)
        Description of additional keyword arguments.

    Returns
    -------
    return_type
        Description of the return value.

    Raises
    ------
    Exception
        Description of the exception.

    Examples
    --------
    >>> function_name(1, 2)
    3       
    """
    pass
```

* [markdownlint](https://github.com/DavidAnson/markdownlint) for markdown documents.
