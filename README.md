# ComFiT: Computational Field Theory Python Package

ComFiT is an open-source Python library for simulating field theories including the Schrödinger equation, damped Gross-Pitaevskii equation, phase-field crystal models, and hydrodynamic models for active nematics and p-adic systems. 
It uses an object-oriented approach to provide functions for system setup, time evolution, data analysis, and visualization. 
The library features the Exponential Time Differencing method for numerical integration, specifically the ETD2RK and ETD4RK schemes, ensuring accurate time evolution in simulations. 

Additionally, ComFiT includes tools for tracking and calculating the density of topological defects in various systems. 
It also enables users to create and export plots and animations for a clearer understanding of the simulated phenomena.

## Features

- **General System Definition**: Utilize the `define_system` superclass to establish general properties of systems.
- **Model-specific Simulations**: Dive deep with subclasses like `PhaseFieldCrystal` for phase-field crystal simulations and `BoseEinsteinCondensate` for Bose-Einstein condensate simulations via the Gross-Pitaevski equation.
- **Extendable Framework**: Easily expand with your own models and simulations.

## Installation

```bash
pip install comfit
```

## Usage

```python
import comfit as cf

# Define a BoseEinsteinCondensate system and simulate
bec= cf.BoseEinsteinCondensate(parameters_here)

# Define a PhaseFieldCrystal system and simulate
pfc = cf.PhaseFieldCrystal(parameters_here)
```

For detailed usage instructions, please refer to our [documentation](/docs/).

## Contributing

We welcome contributions! Please read our [contributing guide](/docs/Contributing.md) to learn about our development process, how to propose bugfixes and improvements, and how to build and test your changes to ComFiT.

## License

ComFiT is licensed under the [MIT License](LICENSE).

## Acknowledgements

We are grateful to [Luiza Angheluta](https://orcid.org/0000-0001-7231-6694) for her steady guidance during our years as Ph.D. doctoral research fellows and introducing us to this field of study, 
[Audun Skaugen](https://orcid.org/0000-0003-0005-786X) for paving the way of the study of these systems in particular and creating the first programs on which this library is built, and 
[Vegard Gjeldvik Jervell](https://orcid.org/0009-0002-2959-0246), for helping us with the technical python parts and becoming full-fledged software developers.  


## Glorified beta-testers

- Harish Pruthviraj Jain
- Milos Joksimovix
