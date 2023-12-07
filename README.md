# ComFiT: Computational Field Theory Python Package

ComFiT is a Python library designed for simulating various condensed matter field theory systems. Developed with expertise from PhD-level research in physics, this package aims to provide a comprehensive toolkit for researchers, students, and enthusiasts in the domain.

## Features

- **General System Definition**: Utilize the `define_system` superclass to establish general properties of systems.
- **Model-specific Simulations**: Dive deep with subclasses like `PFC` for phase-field crystal simulations and `BEC` for Bose-Einstein condensate simulations via the Gross-Pitaevski equation.
- **Extendable Framework**: Easily expand with your own models and simulations.

## Installation

```bash
pip install comfit
```

## Usage

```python
import comfit as cf

# Define a BEC system and simulate
bec= cf.BEC(parameters_here)

# Define a PFC system and simulate
pfc = cf.PFC(parameters_here)
```

For detailed usage instructions, please refer to our [documentation](link_to_docs).

## Setting Up a Development Environment

To ensure a consistent development environment, it's recommended to use a virtual environment. Here's how you can set one up for ComFiT:

### 1. Clone the Repository

First, clone the ComFiT repository to your local machine:

```bash
git clone https://github.com/your-username/ComFiT.git
cd ComFiT
```

Replace `your-username` with the appropriate GitHub username.

### 2. Create a Virtual Environment

You can set up a virtual environment using Python's built-in `venv` module:

```bash
python -m venv venv
```

This creates a `venv` directory in your project root, which will house the virtual environment.

### 3. Activate the Virtual Environment

Depending on your operating system, you can activate the virtual environment as follows:

- **Linux/macOS**:

  ```bash
  source venv/bin/activate
  ```

- **Windows (cmd.exe)**:

  ```bash
  venv\Scripts\activate.bat
  ```

- **Windows (PowerShell)**:

  ```bash
  .\venv\Scripts\Activate.ps1
  ```

Once activated, your terminal or command prompt should indicate that you're inside the `venv` environment.

### 4. Install Dependencies

With the virtual environment activated, you can install the required packages using:

```bash
pip install -r requirements.txt
```

This will install all the dependencies specified in `requirements.txt`.

### 5. Done!

You're now set up with a virtual environment specific to ComFiT. This ensures that any packages or dependencies you install while the environment is activated won't affect your global Python setup.

To deactivate the virtual environment and return to your global Python environment, simply use:

```bash
deactivate
```

## Contributing

We welcome contributions! Please read our [CONTRIBUTING guide](link_to_contributing_guide) to learn about our development process, how to propose bugfixes and improvements, and how to build and test your changes to ComFiT.

## License

ComFiT is licensed under the [MIT License](LICENSE).

## Acknowledgements

We want to express our gratitude to [Institution or Collaborators Names] for their invaluable input and feedback during the development of this package.

## Glorified beta-testers

- Harish Pruthviraj Jain
