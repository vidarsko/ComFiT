# ComFiT: Computational Field Theory Python Package

ComFiT is a Python library designed for simulating various condensed matter field theory systems. Developed with expertise from PhD-level research in physics, this package aims to provide a comprehensive toolkit for researchers, students, and enthusiasts in the domain.

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

For detailed usage instructions, please refer to our [documentation](link_to_docs).

## Contributing

We welcome contributions! Please read our [CONTRIBUTING guide](/docs/Contributing.md) to learn about our development process, how to propose bugfixes and improvements, and how to build and test your changes to ComFiT.

## License

ComFiT is licensed under the [MIT License](LICENSE).

## Acknowledgements

We want to express our gratitude to [Institution or Collaborators Names] for their invaluable input and feedback during the development of this package.

## Glorified beta-testers

- Harish Pruthviraj Jain
