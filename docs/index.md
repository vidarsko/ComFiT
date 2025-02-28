# <img src="img/logo.png" width="25" height="25"> ComFiT documentation

ComFiT ([Github](https://github.com/vidarsko/ComFiT)) is a versatile Python library for simulating field theories, including plotting and animation in an object-oriented manner.
If you use ComFiT in your research, please cite the following paper:

!!! quote ""
    Skogvoll, V., & Rønning, J. (2024). ComFiT: A Python library for computational field theory with topological defects. Journal of Open Source Software, 9(98), 6599. [https://doi.org/10.21105/joss.06599](https://doi.org/10.21105/joss.06599)

Below is a prepromt you can use with a language model to help you get started.

??? abstract "Preprompt for large language model (LLM)"
    ```python
    You are a helpful coding assistant who answers questions to the point.
    
    Gauge the understanding of the user before providing answers.

    Encourage the user to paste error messages if they encounter any.

    Info about ComFiT - Python library for field theories, periodic boundary conditions:

    import comfit as cf (general instance: cfi)

    Class BaseSystem: (instance: bs) (no dynamics)
    Models inheriting from BaseSystem: 
    QuantumMechanics (qm), BoseEinsteinCondensate (bec), NematicLiquidCrystal (nlc), PhaseFieldCrystal (pfc)
    Each model (e.g., qm) directly inherits BaseSystem's attributes, such as dim, dif, and others, accessible directly like, e.g., qm.dif

    Configurable vars:
    dim (1,2, or 3)
    dx
    xmin
    xmax
    xlim ([xmin, xmax])
    xRes
    similar vars for y, z in case bs.dim>1
    dt
    plot_lib ('matplotlib' or 'plotly')

    Other vars: 

    psi (field): primary order parameter (name varies between models)
    psi_f (Fourier transform of psi)
    x (coordinate array)
    xmid
    xmidi (index)
    size_x (xmax-xmin)
    similar vars for y, z in case bs.dim>1
    Res (total)
    dims (xRes if bs.dim=1, [xRes,yRes] if bs.dim=2 etc.)
    rmin = [xmin,ymin,zmin]
    rmid, rmax similar
    volume
    dV
    time (scalar)
    k (list, k[0] wave numbers for x etc.)
    dif (list, dif[i] = 1j*k[i], for differentiation)

    Broadcasting:

    x.shape = (xRes,) if bs.dim=1
    x.shape = (xRes,1) if bs.dim=2, y.shape = (1,yRes)
    similar for x,y,z if bs.dim=3

    Thus, `x+y` is a 2D array of shape `(xRes,yRes)` (no need for meshgrid)

    Functions types:

    calc_-calculates and returns output
    conf_-changes cfi, configures psi and psi_f, returns None
    evolve_-evolves cfi, returns None
    plot_-returns (fig,ax)
    get_-extracts variable 

    Fourier fields denoted `(field)_f`
    Fourier transformation (FT) given by `cfi.fft` and `cfi.ifft`.

    Derivatives using FT, examples:
    
    dxfield = cfi.ifft(bs.dif[0]*field_f)(.real() if field is real)
    Laplacian: cfi.ifft(-bs.calc_k2()*field_f)(.real() if field is real)

    Important functions:

    calc_k2() returns k^2 (for Laplacian)

    Time evolution:
    BaseSystem has no dynamics, but models inheriting BaseSystem have (see below)
    time incremented automatically by `dt` in time evolution loop

    Plotting:
    
    plot_field
    plot_complex_field
    plot_angle_field
    plot_vector_field
    plot_field_in_plane
    plot_complex_field_in_plane
    plot_angle_field_in_plane
    plot_vector_field_in_plane

    Plots (replace `plot_field` under with desired function)

    fig, ax = cfi.plot_field(field, title='title')
    cfi.show(fig) 


    Subplots (if either `number_of_(rows or columns)` is 1, `axs` is list, not 2D array)

    fig, axs = cfi.plot_subplots(2,2)
    cfi.plot_field(field1, fig=fig, ax=axs[0,0])
    cfi.plot_field(field2, fig=fig, ax=axs[0,1]) 
    #etc.

    fig, axs = cfi.plot_subplots(1,2)
    cfi.plot_field(field1, fig=fig, ax=axs[0])
    cfi.plot_field(field2, fig=fig, ax=axs[1]) 
    #etc.

    Animation:

    for n in range(100):
        #Evolve cfi
        fig, ax = cfi.plot_field(field) #replace with appropriate plot function
        cfi.plot_save(n,fig)
    cf.tool_make_animation_gif(n)

    Creating custom model example:

    import comfit as cf
    import numpy as np
    import scipy as sp

    class LandauSystem(cf.BaseSystem):
        def __init__(self,dim, r, **kwargs):
            self.r = r
            super().__init__(dim, **kwargs)
        def calc_omega_f(self):
            return -self.calc_k2() - self.r
        def calc_nonlinear_evolution_function_f(self, field, t):
            return -sp.fft.fftn(field**3) 
        def evolve(self, number_steps):
            omega_f = self.calc_omega_f()
            integrating_factors_f, solver = self.calc_integrating_factors_f_and_solver(omega_f, method='ETD2RK')
            for n in range(number/_steps):
                self.psi, self.psi_f = solver(integrating_factors_f, 
                                            self.calc_nonlinear_evolution_function_f, 
                                            self.psi, self.psi_f)
                self.psi = np.real(self.psi)

    ls = LandauSystem(2, 0.5)
    ls.psi = np.random.rand(ls.xRes,ls.yRes)-0.5
    ls.psi_f = sp.fft.fftn(ls.psi)

    ls.evolve(200)
    fig, ax = ls.plot_field(ls.psi)
    ls.show(fig)

    Models inheriting BaseSystem 

    QuantumMechanics (instance: qm):
    evolve_schrodinger(number_of_steps) evolves qm.psi
    conf_initial_condition_Gaussian(position, width, initial_velocity)
    conf_wavefunction(psi) #sets wavefunction

    BoseEinsteinCondensate (bec)
    evolve_dGPE(number_of_steps) evolves bec.psi
    conf_initial_condition_Thomas_Fermi()
    conf_insert_vortex(charge,position)
    conf_dissipative_frame(interface_width)
    evolve_relax(number_of_steps)
    calc_vortex_nodes()
    plot_nodes(vortex_nodes)

    NematicLiquidCrystal (nlc) Contains 
    evolve_nematic evolves nlc.Q (tensor)
    conf_initial_condition_ordered
    conf_insert_disclination_dipole
    calc_nonlinear_evolution_function_f
    calc_active_force_f
    calc_passive_force_f
    calc_pressure_f
    calc_disclination_density_nematic
    calc_order_and_director
    plot_nodes

    PhaseFieldCrystal (pfc): 
    evolve_PFC
    pfc.psi (real scalar field representing crystalline structures)
    Evolves pfc.psi
    including evolve_PFC. Contains: conf_PFC_from_amplitudes
    calc_PFC_from_amplitudes
    calc_nonlinear_evolution_function_conserved_f
    calc_nonlinear_evolution_function_unconserved_f
    plot_field
    calc_dislocation_nodes
    calc_orientation_field, calc_free_energy
    ```

## Tutorials

The best way to get to know ComFiT is by using it in one of the following tutorials.

###Base System

<div class="grid cards" style="display: flex; flex-wrap: wrap;">
    <a href="https://colab.research.google.com/github/vidarsko/ComFiT/blob/main/tutorial/base_system_basic_framework.ipynb" class="card" style="min-width: 160px; flex: 0 1 calc(20.00% - 10px); margin: 5px;">
        <div> <strong> Basic Framework</strong></div>
        <hr>
        <p>
        <img src="img/index_tutorial_base_system_basic_framework_demo.gif#only-light">
        <img src="img/index_tutorial_base_system_basic_framework_demo-colorinverted.gif#only-dark">
        </p>
        <p style="color: var(--md-default-fg-color)"> Understand the basics of ComFiT, how to calculate derivatives and produce plots and animations in 1, 2 and 3 dimensions. </p>
    </a>
    <a href="https://colab.research.google.com/github/vidarsko/ComFiT/blob/main/tutorial/base_system_make_your_own_model.ipynb" class="card" style="min-width: 160px; flex: 0 1 calc(20.00% - 10px); margin: 5px;">
        <div> <strong>How to make your own model</strong></div>
        <hr>
        <p>
        <img src="img/index_tutorial_base_system_make_your_own_model.gif#only-light">
        <img src="img/index_tutorial_base_system_make_your_own_model-colorinverted.gif#only-dark">
        </p>
        <p style="color: var(--md-default-fg-color)">Learn how to implement, solve and animate your own partial differential equation.</p>
</a>
</div>

### Quantum Mechanics

Modules for learning quantum mechanics with ComFiT.

<div class="grid cards" style="display: flex; flex-wrap: wrap;">
<a href="https://colab.research.google.com/github/vidarsko/ComFiT/blob/main/tutorial/qm_assignment/module_1.ipynb" class="card" style="min-width: 160px; flex: 0 1 calc(20.00% - 10px); margin: 5px;">
    <div> <strong> Module 1 </strong></div>
    <hr>
    <p>
    <img src="img/index_tutorial_qm_assignment_module_1.gif#only-light">
    <img src="img/index_tutorial_qm_assignment_module_1-colorinverted.gif#only-dark">
    </p>
    <p style="color: var(--md-default-fg-color)">
The Schrödinger equation for a single particle, in one and two dimensions.
    </p>
</a>
<a href="https://colab.research.google.com/github/vidarsko/ComFiT/blob/main/tutorial/qm_assignment/module_2.ipynb" class="card" style="min-width: 160px; flex: 0 1 calc(20.00% - 10px); margin: 5px;">
    <div> <strong> Module 2 </strong></div>
    <hr>
    <p>
    <img src="img/index_tutorial_qm_assignment_module_2.png#only-light">
    <img src="img/index_tutorial_qm_assignment_module_2-colorinverted.png#only-dark">
    </p>
    <p style="color: var(--md-default-fg-color)">
    Operators and expectation values.
    </p>
</a>
<a href="https://colab.research.google.com/github/vidarsko/ComFiT/blob/main/tutorial/qm_assignment/module_3.ipynb" class="card" style="min-width: 160px; flex: 0 1 calc(20.00% - 10px); margin: 5px;">
    <div> <strong> Module 3 </strong></div>
    <hr>
    <p>
    <img src="img/index_tutorial_qm_harmonic_oscillator.gif#only-light">
    <img src="img/index_tutorial_qm_harmonic_oscillator-colorinverted.gif#only-dark">
    </p>
    <p style="color: var(--md-default-fg-color)">
    The Quantum Harmonic Oscillator and her eigenstates.
    </p>
</a>
</div>

<!-- 
<div class="grid cards" style="display: flex; flex-wrap: wrap;">
    <a href="https://colab.research.google.com/github/vidarsko/ComFiT/blob/main/tutorial/quantum_mechanics_1D_wave_packet.ipynb" class="card" style="min-width: 160px; flex: 0 1 calc(20.00% - 10px); margin: 5px;">
        <div> <strong> 1D wave packets </strong></div>
        <hr>
        <p>
        <img src="img/index_1D_Quantum_Mechanics_Gaussian.gif#only-light">
        <img src="img/index_1D_Quantum_Mechanics_Gaussian-colorinverted.gif#only-dark">
        </p>
        <p style="color: var(--md-default-fg-color)">
        Understand the basics of the Quantum Mechanics model.
        </p>
    </a>
    <a href="https://colab.research.google.com/github/vidarsko/ComFiT/blob/main/tutorial/quantum_mechanics_2D_wave_packet.ipynb" class="card" style="min-width: 160px; flex: 0 1 calc(20.00% - 10px); margin: 5px;">
        <div> <strong> 2D wave packets </strong></div>
        <hr>
        <p>
        <img src="img/index_tutorial_qm_2D_wave_packet.gif#only-light">
        <img src="img/index_tutorial_qm_2D_wave_packet-colorinverted.gif#only-dark">
        </p>
        <p style="color: var(--md-default-fg-color)">
        Understand how to plot a quantum mechanical system in 2 dimensions.
        </p>
</a>
<a href="https://colab.research.google.com/github/vidarsko/ComFiT/blob/main/tutorial/quantum_mechanics_3D_wave_packet.ipynb" class="card" style="min-width: 160px; flex: 0 1 calc(20.00% - 10px); margin: 5px;">
    <div> <strong> 3D wave packets </strong></div>
    <hr>
    <p>
    <img src="img/index_tutorial_qm_3D_wave_packet.gif#only-light">
    <img src="img/index_tutorial_qm_3D_wave_packet-colorinverted.gif#only-dark">
    </p>
    <p style="color: var(--md-default-fg-color)">
    Understand how to plot a quantum mechanical system in 3 dimensions.
    </p>
</a>
</a>
<a href="https://colab.research.google.com/github/vidarsko/ComFiT/blob/main/tutorial/quantum_mechanics_harmonic_oscillator.ipynb" class="card" style="min-width: 160px; flex: 0 1 calc(20.00% - 10px); margin: 5px;">
    <div> <strong> The harmonic oscillator </strong></div>
    <hr>
    <p>
    <img src="img/index_tutorial_qm_harmonic_oscillator.gif#only-light">
    <img src="img/index_tutorial_qm_harmonic_oscillator-colorinverted.gif#only-dark">
    </p>
    <p style="color: var(--md-default-fg-color)">
    Investigate the quantum mechanical harmonic oscillator.
    </p>
</a>
<a href="https://colab.research.google.com/github/vidarsko/ComFiT/blob/main/tutorial/quantum_mechanics_the_hydrogen_atom.ipynb" class="card" style="min-width: 160px; flex: 0 1 calc(20.00% - 10px); margin: 5px;">
    <div> <strong> The hydrogen atom </strong></div>
    <hr>
    <p>
    <img src="img/index_tutorial_quantum_mechanics_hydrogen.gif#only-light">
    <img src="img/index_tutorial_quantum_mechanics_hydrogen-colorinverted.gif#only-dark">
    </p>
    <p style="color: var(--md-default-fg-color)">
    Get to know the hydrogen atom.
    </p>
</a>
</div> -->

### Bose-Einstein Condensates

<div class="grid cards" style="display: flex; flex-wrap: wrap;">
<a href="https://colab.research.google.com/github/vidarsko/ComFiT/blob/main/tutorial/bose_einstein_condensate_basic_framework.ipynb" class="card" style="min-width: 160px; flex: 0 1 calc(20.00% - 10px); margin: 5px;">
    <div> <strong> Basic Framework </strong></div>
    <hr>
    <p>
    <img src="img/index_tutorial_bec_basic_framework-colorinverted.gif#only-dark">
    <img src="img/index_tutorial_bec_basic_framework.gif#only-light">
    </p>
    <p style="color: var(--md-default-fg-color)">
    Understand the basics of the Bose Einstein Condensate model.
    </p>
</a>
<a href="https://colab.research.google.com/github/vidarsko/ComFiT/blob/main/tutorial/bose_einstein_condensate_time_dependenent_potentials.ipynb" class="card" style="min-width: 160px; flex: 0 1 calc(20.00% - 10px); margin: 5px;">
    <div> <strong>Time-dependent potentials</strong></div>
    <hr>
    <p>
    <img src="img/index_tutorial_bec_time_dependent_potentials-colorinverted.gif#only-dark">
    <img src="img/index_tutorial_bec_time_dependent_potentials.gif#only-light">
    </p>
    <p style="color: var(--md-default-fg-color)">
    Learn how to create time-dependent potentials to stir the Bose Einstein Condensate model.
    </p>
</a>
<a href="https://colab.research.google.com/github/vidarsko/ComFiT/blob/main/tutorial/bose_einstein_condensate_comoving_frame_and_defect_tracking.ipynb" class="card" style="min-width: 160px; flex: 0 1 calc(20.00% - 10px); margin: 5px;">
    <div> <strong>Comoving frame and defect tracking</strong></div>
    <hr>
    <p>
    <img src="img/index_tutorial_bec_comoving_frame_defect_tracking-colorinverted.gif#only-dark">
    <img src="img/index_tutorial_bec_comoving_frame_defect_tracking.gif#only-light">
    </p>
    <p style="color: var(--md-default-fg-color)">
    Learn how to track defects and study defects made by an obstacle.
    </p>
</a>
<a href="https://colab.research.google.com/github/vidarsko/ComFiT/blob/main/tutorial/bose_einstein_condensate_3D_comoving_frame.ipynb" class="card" style="min-width: 160px; flex: 0 1 calc(20.00% - 10px); margin: 5px;">
    <div> <strong>3D and comoving frame </strong></div>
    <hr>
    <p>
    <img src="img/index_tutorial_bec_3D_comoving_frame-colorinverted.gif#only-dark">
    <img src="img/index_tutorial_bec_3D_comoving_frame.gif#only-light">
    </p>
    <p style="color: var(--md-default-fg-color)">
    Learn how to use the Bose Einstein Condensate model in 3 dimensions and in a comoving frame.
    </p>
</a>
</div>

### Nematic Liquid Crystal

<div class="grid cards" style="display: flex; flex-wrap: wrap;">
<a href="https://colab.research.google.com/github/vidarsko/ComFiT/blob/main/tutorial/nematic_liquid_crystal_2D_active_nematic.ipynb" class="card" style="min-width: 160px; flex: 0 1 calc(20.00% - 10px); margin: 5px;">
    <div> <strong> 2D active nematic </strong></div>
    <hr>
    <p>
    <img src="img/index_tutorial_nematic_liquid_crystal_2D_active_nematic.gif#only-light">
    <img src="img/index_tutorial_nematic_liquid_crystal_2D_active_nematic-colorinverted.gif#only-dark">
    </p>
    <p style="color: var(--md-default-fg-color)">
    Simulate an active nematic liquid crystal in 2 dimensions.
    </p>
</a>
</div>

### Phase-field crystal

<div class="grid cards" style="display: flex; flex-wrap: wrap;">
<a href="https://colab.research.google.com/github/vidarsko/ComFiT/blob/main/tutorial/phase_field_crystal_basic_framework.ipynb" class="card" style="min-width: 160px; flex: 0 1 calc(20.00% - 10px); margin: 5px;">
    <div> <strong> Basic framework </strong></div>
    <hr>
    <p>
    <img src="img/index_tutorial_pfc_basic_framework.gif#only-light">
    <img src="img/index_tutorial_pfc_basic_framework-colorinverted.gif#only-dark">
    </p>
    <p style="color: var(--md-default-fg-color)">
    Get to know the basic of the PFC framework, including how to insert dislocations, plot them, and evolve the PFC.
    </p>
<a href="https://colab.research.google.com/github/vidarsko/ComFiT/blob/main/tutorial/phase_field_crystal_stresses_and_strains.ipynb" class="card" style="min-width: 160px; flex: 0 1 calc(20.00% - 10px); margin: 5px;">
    <div> <strong> Stresses and strains </strong></div>
    <hr>
    <p>
    <img src="img/index_tutorial_pfc_stresses_and_strains.gif#only-light">
    <img src="img/index_tutorial_pfc_stresses_and_strains-colorinverted.gif#only-dark">
    </p>
    <p style="color: var(--md-default-fg-color)">
    Learn how to calculate stresses and strains in the PFC model.
    </p>
</a>
<a href="https://colab.research.google.com/github/vidarsko/ComFiT/blob/main/tutorial/phase_field_crystal_polycrystalline_systems.ipynb" class="card" style="min-width: 160px; flex: 0 1 calc(20.00% - 10px); margin: 5px;">
    <div> <strong> Polycrystalline systems </strong></div>
    <hr>
    <p>
    <img src="img/index_tutorial_pfc_polycrystals.gif#only-light">
    <img src="img/index_tutorial_pfc_polycrystals-colorinverted.gif#only-dark">
    </p>
    <p style="color: var(--md-default-fg-color)">
    Create polycrystalline systems with the PFC model and evolve it according to different dynamics.
    </p>
</a>
</div>


For the time being, ComFiT is limited to periodic boundary conditions, but this may change in the future.

## Installation

Comfit can be installed from the Python Package Index (PyPI), a repository of software for the Python programming language, by executing the command

```bash
pip install comfit
```

pip install comfit in your terminal or command prompt.

## Virtual environnement

Using a virtual environnement when using ComFiT is highly encouraged for because even though we try to write robust code, it is still a library under development, so previously written simulations may break. By keeping your simulations together with the specific version of ComFiT, you make sure that your simulations will not break due to coming updates.

To create a virtual environnement, run the following command in your terminal after having navigated to the root folder of your exploration project

```bash
Python -m venv myvenv
```

This will create the folder `myvenv` which will contain the local installation of Python and associated packages.
To activate the virtual environnement, simply run

```bash
.\venv\Scripts\activate
```

from the terminal.
Afterwards, you may install ComFiT using PyPi.
If your folder is part of a github repository, it is recommended to remove the virtual environment from the git project by adding `venv/` to your `.gitignore` file.

## Contributing

We welcome contributions.
Whether you're fixing a bug, adding a new feature, or improving our documentation, your support helps us make the package more robust and versatile.
Contributions can take many forms, from fixing minor bugs to implementing complex new features. 
Below are the ways you can contribute:

### Bug Fixes

Did you identify a bug? Here's how to proceed:

1. **Fork the repository**: Start by forking the ComFiT GitHub repository.
2. **Create a branch**: Make a new branch on your fork dedicated to the bug fix.
3. **Fix the bug**: Make the necessary changes to resolve the bug.
4. **Run tests**: Ensure all existing tests pass with your changes. Add new tests if necessary to cover the bug fix.
5. **Submit a Pull Request (PR)**: Create a PR against the main ComFiT repository. Clearly describe the bug and how your changes fix it.

### Reporting Issues

Encountered an issue or have a suggestion? Please follow these steps:

1. **Check existing issues**: Before creating a new issue, please check existing issues to avoid duplicates.
2. **Create a new issue**: If your issue is unique, open a new issue on GitHub. Provide a detailed description, including steps to reproduce the issue if applicable.

### Feature Requests

Got an idea for a new feature or enhancement? We'd love to hear it! Please raise a discussion, or an issue as outlined above, detailing your idea and its potential benefits to ComFiT.

### Adding Your Own Model

If you're interested in adding your own model to ComFiT, we welcome your contribution! Your model should adhere to the following guidelines:

1. **Well-documented**: Include detailed documentation explaining your model's theory, implementation, and usage.
2. **Thoroughly tested**: Write comprehensive tests covering the functionality of your model.
3. **Follow ComFiT structure**: Ensure your model integrates seamlessly with the existing ComFiT framework.
4. **Tutorial**: Consider adding a tutorial in the form of a Jupyter notebook, demonstrating how to use your model. Link to the tutorial in your contribution.

For detailed instructions on implementing your own PDE model with ComFiT, refer to our [tutorial for creating your own model](https://colab.research.google.com/github/vidarsko/ComFiT/blob/main/tutorial/base_system_make_your_own_model.ipynb).

### Documentation Improvements

Good documentation is key to a project's usability and its community's growth. 
If you see areas for improvement or want to add documentation for undocumented features, your contributions are greatly appreciated.
