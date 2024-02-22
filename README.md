# ComFiT: Computational Field Theory Python Package

ComFiT is an open-source Python library for simulating field theories including the Schrödinger equation, damped Gross-Pitaevskii equation, phase-field crystal models, and hydrodynamic models for active nematics and p-adic systems.
It uses an object-oriented approach to provide functions for system setup, time evolution, data analysis, and visualization.
The library features the Exponential Time Differencing method for numerical integration, specifically the ETD2RK and ETD4RK schemes, ensuring accurate time evolution in simulations.

Additionally, ComFiT includes tools for tracking and calculating the density of topological defects in various systems.
It also enables users to create and export plots and animations for a clearer understanding of the simulated phenomena.

Below, we see the simulation of an initial Guassian wave function in a harmonic oscillator trap under the time evolution imposed by the Schrödinger equation.

![Oscillating 1D Gaussian wave packet in a harmonic oscillator](docs/images/index_1D_Quantum_Mechanics_Gaussian.gif)

*Figure: Oscillating 1D Gaussian wave packet in a harmonic oscillator* A 1D simulation of the Schrödinger equation with a Gaussian initial condition (at rest) in a Harmonic oscillator potential.

![Oscillating 3D Gaussian wave packet in a harmonic oscillator](docs/images/index_3D_Quantum_Mechanics_Gaussian.gif)

*Figure: Oscillating 3D Gaussian wave packet in a harmonic oscillator* A 3D simulation of the Schrödinger equation with a Gaussian initial condition (at rest) in a Harmonic oscillator potential.
The amplitude of the wavefunciton is shown as an isosurface at $50\%$ its maximal value and the argument by the angle color scheme of the figure above.

## Documentation

For detailed usage instructions, tutorials and instructions on how to contribute, please refer to our [documentation](https://vidarsko.github.io/ComFiT/).

## Installation

Comfit can be installed from the Python Package Index (PyPI), a repository of software for the Python programming language, by executing the command

```bash
pip install comfit
```

pip install comfit in your terminal or command prompt.
Using a virtual environment is highly recommended as new comfit is an open source project subject to likely frequent updates and overall bug fixes.

## License

ComFiT is licensed under the [MIT License](LICENSE).

## Acknowledgements

We are grateful to [Luiza Angheluta](https://orcid.org/0000-0001-7231-6694) for her steady guidance during our years as Ph.D. doctoral research fellows and introducing us to this field of study,
[Audun Skaugen](https://orcid.org/0000-0003-0005-786X) for paving the way of the study of these systems in particular and creating the first programs on which this library is built, and
[Vegard Gjeldvik Jervell](https://orcid.org/0009-0002-2959-0246), for helping us with the technical python parts and becoming full-fledged software developers.
