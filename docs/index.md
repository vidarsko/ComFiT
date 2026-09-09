# <img src="img/logo.png" width="25" height="25" alt="ComFiT logo"> ComFiT documentation

ComFiT ([Github](https://github.com/vidarsko/ComFiT)) is a versatile Python library for simulating field theories, including plotting and animation in an object-oriented manner.
If you use ComFiT in your research, please cite the following paper:

!!! quote ""
    Skogvoll, V., & Rønning, J. (2024). ComFiT: A Python library for computational field theory with topological defects. Journal of Open Source Software, 9(98), 6599. [https://doi.org/10.21105/joss.06599](https://doi.org/10.21105/joss.06599)

Below is a preprompt you can paste into a language model's chat window (before your own question) to help it give better answers about ComFiT.

<!-- markdownlint-disable MD046 -->
??? abstract "Preprompt for large language model (LLM)"
    ```python
    You are a helpful coding assistant answering questions about ComFiT, a Python library for
    simulating field theories with periodic boundary conditions.

    Answer to the point. Gauge the user's level of understanding before going into depth.
    You may not always be correct about ComFiT's API - if unsure, say so, and ask the user to
    paste the exact error message/traceback if they hit one.

    ## Package structure

    import comfit as cf (general instance name convention: cfi)

    BaseSystem (instance: bs): the base class, holds the grid/Fourier machinery, no dynamics of
    its own. All models below inherit from it, so any BaseSystem attribute or method (e.g. bs.dif,
    bs.fft) is also available directly on qm/bec/nlc/pfc instances.

    Models inheriting from BaseSystem:
    QuantumMechanics (qm), BoseEinsteinCondensate (bec), NematicLiquidCrystal (nlc),
    PhaseFieldCrystal (pfc, abstract - instantiate one of PhaseFieldCrystal1DPeriodic,
    PhaseFieldCrystal2DTriangular, PhaseFieldCrystal2DSquare,
    PhaseFieldCrystal3DBodyCenteredCubic, PhaseFieldCrystal3DFaceCenteredCubic,
    PhaseFieldCrystal3DSimpleCubic)

    ## Configuration (constructor kwargs, also readable as attributes afterwards)

    dim (1, 2, or 3)
    dx, xmin, xmax, xlim ([xmin, xmax]), xRes - and the equivalent y*/z* variants when bs.dim > 1
    dt
    plot_lib ('matplotlib' or 'plotly', default 'plotly')

    ## Other attributes

    psi (field): primary order parameter (name/meaning varies by model, e.g. bec.psi is the
    condensate wavefunction, pfc.psi is the crystal density field, nlc.Q is a tensor field instead)
    psi_f: Fourier transform of psi
    x: coordinate array (and y, z when bs.dim > 1)
    xmid, xmidi: midpoint coordinate and its index (and y/z equivalents)
    size_x: xmax - xmin (and y/z equivalents)
    Res: total number of grid points
    dims: xRes if bs.dim==1, [xRes,yRes] if bs.dim==2, etc.
    rmin = [xmin,ymin,zmin], rmid and rmax similar
    volume, dV: cell volume element
    time: current simulation time (scalar), advanced by dt each evolve_* step
    k: list of wavenumber arrays, k[0] for x etc.
    dif: list of spectral derivative operators, dif[i] = 1j*k[i]

    ## Broadcasting

    x.shape = (xRes,) if bs.dim==1
    x.shape = (xRes,1) if bs.dim==2, y.shape = (1,yRes)
    similarly for x,y,z if bs.dim==3, so `x+y` broadcasts to shape (xRes,yRes) with no meshgrid
    needed

    ## Method name prefixes

    calc_ - computes and returns a value, no side effects
    conf_ - configures/mutates the instance (e.g. sets psi and psi_f), returns None
    evolve_ - advances the instance in time, returns None
    plot_ - returns (fig, ax)
    get_ - extracts/reads a derived variable

    Fourier-space fields are named `<field>_f`. Transform with `cfi.fft` / `cfi.ifft`.

    Spectral derivatives, examples:

    dxfield = cfi.ifft(cfi.dif[0] * field_f)  # .real if field is a real-valued field
    laplacian = cfi.ifft(-cfi.calc_k2() * field_f)  # calc_k2() returns k^2

    ## Time evolution

    BaseSystem has no dynamics; each model implements its own evolve_* method(s) (see below).
    Calling evolve_*(number_of_steps) advances psi (or the model's field) and increments
    cfi.time automatically.

    ## Plotting

    plot_field, plot_complex_field, plot_angle_field, plot_vector_field, and the *_in_plane
    variants of each (plot_field_in_plane, plot_complex_field_in_plane, etc., for slicing 3D
    fields)

    fig, ax = cfi.plot_field(field, title='title')
    cfi.show(fig)

    Subplots (axs is a flat list, not a 2D array, if either grid dimension is 1):

    fig, axs = cfi.plot_subplots(2, 2)
    cfi.plot_field(field1, fig=fig, ax=axs[0, 0])
    cfi.plot_field(field2, fig=fig, ax=axs[0, 1])
    # etc.

    fig, axs = cfi.plot_subplots(1, 2)
    cfi.plot_field(field1, fig=fig, ax=axs[0])
    cfi.plot_field(field2, fig=fig, ax=axs[1])
    # etc.

    Animation:

    number_of_frames = 100
    for n in range(number_of_frames):
        # evolve cfi here
        fig, ax = cfi.plot_field(field)  # replace with the appropriate plot function
        cfi.plot_save(fig, n)
    cf.tool_make_animation_gif(number_of_frames - 1)

    ## Creating a custom model

    import comfit as cf
    import numpy as np
    import scipy as sp

    class LandauSystem(cf.BaseSystem):
        def __init__(self, dim, r, **kwargs):
            self.r = r
            super().__init__(dim, **kwargs)
        def calc_omega_f(self):
            return -self.calc_k2() - self.r
        def calc_nonlinear_evolution_function_f(self, field, t):
            return -sp.fft.fftn(field**3)
        def evolve(self, number_of_steps):
            omega_f = self.calc_omega_f()
            integrating_factors_f, solver = self.calc_integrating_factors_f_and_solver(omega_f, method='ETD2RK')
            for n in range(number_of_steps):
                self.psi, self.psi_f = solver(integrating_factors_f,
                                            self.calc_nonlinear_evolution_function_f,
                                            self.psi, self.psi_f)
                self.psi = np.real(self.psi)

    ls = LandauSystem(2, 0.5)
    ls.psi = np.random.rand(ls.xRes, ls.yRes) - 0.5
    ls.psi_f = sp.fft.fftn(ls.psi)

    ls.evolve(200)
    fig, ax = ls.plot_field(ls.psi)
    ls.show(fig)

    ## Models inheriting BaseSystem

    QuantumMechanics (instance: qm) - qm.psi is the wavefunction:
    evolve_schrodinger(number_of_steps) evolves qm.psi
    conf_initial_condition_gaussian(position, width, initial_velocity)
    conf_wavefunction(psi) # sets wavefunction directly

    BoseEinsteinCondensate (bec) - bec.psi is the condensate wavefunction:
    evolve_dGPE(number_of_steps) evolves bec.psi
    evolve_relax(number_of_steps) relaxes bec.psi towards a stationary state
    conf_initial_condition_thomas_fermi()
    conf_insert_vortex(charge, position)
    conf_dissipative_frame(interface_width)
    calc_vortex_nodes() returns detected vortex node positions/charges
    plot_nodes(vortex_nodes)

    NematicLiquidCrystal (nlc) - nlc.Q is the tensor order parameter:
    evolve_nematic(number_of_steps) evolves nlc.Q
    conf_initial_condition_ordered
    conf_insert_disclination_dipole
    calc_active_force_f, calc_passive_force_f, calc_pressure_f
    calc_disclination_density
    calc_order_and_director
    calc_disclination_nodes() returns detected disclination node positions/charges
    plot_nodes

    PhaseFieldCrystal (pfc) - pfc.psi is the (real, scalar) crystal density field:
    evolve_PFC(number_of_steps) evolves pfc.psi
    conf_PFC_from_amplitudes / calc_PFC_from_amplitudes
    calc_nonlinear_evolution_function_conserved_f, calc_nonlinear_evolution_function_unconserved_f
    calc_dislocation_nodes() returns detected dislocation node positions/Burgers vectors
    calc_orientation_field, calc_free_energy
    plot_field
    ```
<!-- markdownlint-enable MD046 -->

See the ComFiT Library Reference below for a complete list of class methods and their usage.

<div class="grid cards" style="display: flex; flex-wrap: wrap;">
    <a href="https://comfitlib.com/library_reference/" class="card" style="min-width: 160px; flex: 0 1 calc(100.00% - 10px); margin: 5px;">
        <div style="text-align: center;">
            <strong> ComFiT Library Reference </strong>
        </div>
    </a>
</div>

## Tutorials

The best way to get to know ComFiT is by using it in one of the following tutorials.

### Base System

<div class="grid cards" style="display: flex; flex-wrap: wrap;">
    <a href="https://colab.research.google.com/github/vidarsko/ComFiT/blob/main/tutorial/base_system_basic_framework.ipynb" class="card" style="min-width: 160px; flex: 0 1 calc(20.00% - 10px); margin: 5px;">
        <div> <strong> Basic Framework</strong></div>
        <hr>
        <p>
        <img src="img/index_tutorial_base_system_basic_framework_demo.gif#only-light" alt="Computing derivatives and animating fields with BaseSystem">
        <img src="img/index_tutorial_base_system_basic_framework_demo-colorinverted.gif#only-dark" alt="Computing derivatives and animating fields with BaseSystem">
        </p>
        <p style="color: var(--md-default-fg-color)"> Understand the basics of ComFiT, how to calculate derivatives and produce plots and animations in 1, 2 and 3 dimensions. </p>
    </a>
    <a href="https://colab.research.google.com/github/vidarsko/ComFiT/blob/main/tutorial/base_system_make_your_own_model.ipynb" class="card" style="min-width: 160px; flex: 0 1 calc(20.00% - 10px); margin: 5px;">
        <div> <strong>How to make your own model</strong></div>
        <hr>
        <p>
        <img src="img/index_tutorial_base_system_make_your_own_model.gif#only-light" alt="Implementing and animating a custom partial differential equation model">
        <img src="img/index_tutorial_base_system_make_your_own_model-colorinverted.gif#only-dark" alt="Implementing and animating a custom partial differential equation model">
        </p>
        <p style="color: var(--md-default-fg-color)">Learn how to implement, solve and animate your own partial differential equation.</p>
</a>
</div>

### Quantum Mechanics

Modules for learning quantum mechanics with ComFiT.
Author: [Carl Fredrik Nordbø Knutsen](https://www.mn.uio.no/fysikk/?vrtx=person-view&uid=cfknutse).

<div class="grid cards" style="display: flex; flex-wrap: wrap;">
<a href="https://colab.research.google.com/github/vidarsko/ComFiT/blob/main/tutorial/qm_assignment/module_1.ipynb" class="card" style="min-width: 160px; flex: 0 1 calc(20.00% - 10px); margin: 5px;">
    <div> <strong> Module 1 </strong></div>
    <hr>
    <p>
    <img src="img/index_tutorial_qm_assignment_module_1.gif#only-light" alt="Schrödinger equation simulation for a single particle in one and two dimensions">
    <img src="img/index_tutorial_qm_assignment_module_1-colorinverted.gif#only-dark" alt="Schrödinger equation simulation for a single particle in one and two dimensions">
    </p>
    <p style="color: var(--md-default-fg-color)">
The Schrödinger equation for a single particle, in one and two dimensions.
    </p>
</a>
<a href="https://colab.research.google.com/github/vidarsko/ComFiT/blob/main/tutorial/qm_assignment/module_2.ipynb" class="card" style="min-width: 160px; flex: 0 1 calc(20.00% - 10px); margin: 5px;">
    <div> <strong> Module 2 </strong></div>
    <hr>
    <p>
    <img src="img/index_tutorial_qm_assignment_module_2.png#only-light" alt="Quantum operators and expectation values">
    <img src="img/index_tutorial_qm_assignment_module_2-colorinverted.png#only-dark" alt="Quantum operators and expectation values">
    </p>
    <p style="color: var(--md-default-fg-color)">
    Operators and expectation values.
    </p>
</a>
<a href="https://colab.research.google.com/github/vidarsko/ComFiT/blob/main/tutorial/qm_assignment/module_3.ipynb" class="card" style="min-width: 160px; flex: 0 1 calc(20.00% - 10px); margin: 5px;">
    <div> <strong> Module 3 </strong></div>
    <hr>
    <p>
    <img src="img/index_tutorial_qm_harmonic_oscillator.gif#only-light" alt="Quantum harmonic oscillator eigenstates">
    <img src="img/index_tutorial_qm_harmonic_oscillator-colorinverted.gif#only-dark" alt="Quantum harmonic oscillator eigenstates">
    </p>
    <p style="color: var(--md-default-fg-color)">
    The Quantum Harmonic Oscillator and her eigenstates.
    </p>
</a>
<a href="https://colab.research.google.com/github/vidarsko/ComFiT/blob/main/tutorial/qm_assignment/module_4.ipynb" class="card" style="min-width: 160px; flex: 0 1 calc(20.00% - 10px); margin: 5px;">
    <div> <strong> Module 4 </strong></div>
    <hr>
    <p>
    <img src="img/quantum_mechanics_barrier_reflection.gif#only-light" alt="Quantum tunneling through a potential barrier">
    <img src="img/quantum_mechanics_barrier_reflection-colorinverted.gif#only-dark" alt="Quantum tunneling through a potential barrier">
    </p>
    <p style="color: var(--md-default-fg-color)">
    Quantum tunneling.
    </p>
</a>
</div>

General tutorials.

<div class="grid cards" style="display: flex; flex-wrap: wrap;">
    <a href="https://colab.research.google.com/github/vidarsko/ComFiT/blob/main/tutorial/quantum_mechanics_1D_wave_packet.ipynb" class="card" style="min-width: 160px; flex: 0 1 calc(20.00% - 10px); margin: 5px;">
        <div> <strong> 1D wave packets </strong></div>
        <hr>
        <p>
        <img src="img/index_1D_Quantum_Mechanics_Gaussian.gif#only-light" alt="1D Gaussian wave packet simulation">
        <img src="img/index_1D_Quantum_Mechanics_Gaussian-colorinverted.gif#only-dark" alt="1D Gaussian wave packet simulation">
        </p>
        <p style="color: var(--md-default-fg-color)">
        Understand the basics of the Quantum Mechanics model.
        </p>
    </a>
    <a href="https://colab.research.google.com/github/vidarsko/ComFiT/blob/main/tutorial/quantum_mechanics_2D_wave_packet.ipynb" class="card" style="min-width: 160px; flex: 0 1 calc(20.00% - 10px); margin: 5px;">
        <div> <strong> 2D wave packets </strong></div>
        <hr>
        <p>
        <img src="img/index_tutorial_qm_2D_wave_packet.gif#only-light" alt="2D quantum wave packet simulation">
        <img src="img/index_tutorial_qm_2D_wave_packet-colorinverted.gif#only-dark" alt="2D quantum wave packet simulation">
        </p>
        <p style="color: var(--md-default-fg-color)">
        Understand how to plot a quantum mechanical system in 2 dimensions.
        </p>
</a>
<a href="https://colab.research.google.com/github/vidarsko/ComFiT/blob/main/tutorial/quantum_mechanics_3D_wave_packet.ipynb" class="card" style="min-width: 160px; flex: 0 1 calc(20.00% - 10px); margin: 5px;">
    <div> <strong> 3D wave packets </strong></div>
    <hr>
    <p>
    <img src="img/index_tutorial_qm_3D_wave_packet.gif#only-light" alt="3D quantum wave packet simulation">
    <img src="img/index_tutorial_qm_3D_wave_packet-colorinverted.gif#only-dark" alt="3D quantum wave packet simulation">
    </p>
    <p style="color: var(--md-default-fg-color)">
    Understand how to plot a quantum mechanical system in 3 dimensions.
    </p>
</a>
<!-- </a>
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
</a> -->
<a href="https://colab.research.google.com/github/vidarsko/ComFiT/blob/main/tutorial/quantum_mechanics_the_hydrogen_atom.ipynb" class="card" style="min-width: 160px; flex: 0 1 calc(20.00% - 10px); margin: 5px;">
    <div> <strong> The hydrogen atom </strong></div>
    <hr>
    <p>
    <img src="img/index_tutorial_quantum_mechanics_hydrogen.gif#only-light" alt="Hydrogen atom wavefunction simulation">
    <img src="img/index_tutorial_quantum_mechanics_hydrogen-colorinverted.gif#only-dark" alt="Hydrogen atom wavefunction simulation">
    </p>
    <p style="color: var(--md-default-fg-color)">
    Get to know the hydrogen atom.
    </p>
</a>
</div>

### Bose-Einstein Condensates

<div class="grid cards" style="display: flex; flex-wrap: wrap;">
<a href="https://colab.research.google.com/github/vidarsko/ComFiT/blob/main/tutorial/bose_einstein_condensate_basic_framework.ipynb" class="card" style="min-width: 160px; flex: 0 1 calc(20.00% - 10px); margin: 5px;">
    <div> <strong> Basic Framework </strong></div>
    <hr>
    <p>
    <img src="img/index_tutorial_bec_basic_framework-colorinverted.gif#only-dark" alt="Basic Bose-Einstein condensate simulation framework">
    <img src="img/index_tutorial_bec_basic_framework.gif#only-light" alt="Basic Bose-Einstein condensate simulation framework">
    </p>
    <p style="color: var(--md-default-fg-color)">
    Understand the basics of the Bose Einstein Condensate model.
    </p>
</a>
<a href="https://colab.research.google.com/github/vidarsko/ComFiT/blob/main/tutorial/bose_einstein_condensate_time_dependenent_potentials.ipynb" class="card" style="min-width: 160px; flex: 0 1 calc(20.00% - 10px); margin: 5px;">
    <div> <strong>Time-dependent potentials</strong></div>
    <hr>
    <p>
    <img src="img/index_tutorial_bec_time_dependent_potentials-colorinverted.gif#only-dark" alt="Stirring a Bose-Einstein condensate with a time-dependent potential">
    <img src="img/index_tutorial_bec_time_dependent_potentials.gif#only-light" alt="Stirring a Bose-Einstein condensate with a time-dependent potential">
    </p>
    <p style="color: var(--md-default-fg-color)">
    Learn how to create time-dependent potentials to stir the Bose Einstein Condensate model.
    </p>
</a>
<a href="https://colab.research.google.com/github/vidarsko/ComFiT/blob/main/tutorial/bose_einstein_condensate_comoving_frame_and_defect_tracking.ipynb" class="card" style="min-width: 160px; flex: 0 1 calc(20.00% - 10px); margin: 5px;">
    <div> <strong>Comoving frame and defect tracking</strong></div>
    <hr>
    <p>
    <img src="img/index_tutorial_bec_comoving_frame_defect_tracking-colorinverted.gif#only-dark" alt="Tracking vortex defects in a Bose-Einstein condensate in a comoving frame">
    <img src="img/index_tutorial_bec_comoving_frame_defect_tracking.gif#only-light" alt="Tracking vortex defects in a Bose-Einstein condensate in a comoving frame">
    </p>
    <p style="color: var(--md-default-fg-color)">
    Learn how to track defects and study defects made by an obstacle.
    </p>
</a>
<a href="https://colab.research.google.com/github/vidarsko/ComFiT/blob/main/tutorial/bose_einstein_condensate_3D_comoving_frame.ipynb" class="card" style="min-width: 160px; flex: 0 1 calc(20.00% - 10px); margin: 5px;">
    <div> <strong>3D and comoving frame </strong></div>
    <hr>
    <p>
    <img src="img/index_tutorial_bec_3D_comoving_frame-colorinverted.gif#only-dark" alt="3D Bose-Einstein condensate simulation in a comoving frame">
    <img src="img/index_tutorial_bec_3D_comoving_frame.gif#only-light" alt="3D Bose-Einstein condensate simulation in a comoving frame">
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
    <img src="img/index_tutorial_nematic_liquid_crystal_2D_active_nematic.gif#only-light" alt="2D active nematic liquid crystal simulation">
    <img src="img/index_tutorial_nematic_liquid_crystal_2D_active_nematic-colorinverted.gif#only-dark" alt="2D active nematic liquid crystal simulation">
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
    <img src="img/index_tutorial_pfc_basic_framework.gif#only-light" alt="Basic phase field crystal framework with dislocations">
    <img src="img/index_tutorial_pfc_basic_framework-colorinverted.gif#only-dark" alt="Basic phase field crystal framework with dislocations">
    </p>
    <p style="color: var(--md-default-fg-color)">
    Get to know the basic of the PFC framework, including how to insert dislocations, plot them, and evolve the PFC.
    </p>
<a href="https://colab.research.google.com/github/vidarsko/ComFiT/blob/main/tutorial/phase_field_crystal_stresses_and_strains.ipynb" class="card" style="min-width: 160px; flex: 0 1 calc(20.00% - 10px); margin: 5px;">
    <div> <strong> Stresses and strains </strong></div>
    <hr>
    <p>
    <img src="img/index_tutorial_pfc_stresses_and_strains.gif#only-light" alt="Calculating stresses and strains in a phase field crystal model">
    <img src="img/index_tutorial_pfc_stresses_and_strains-colorinverted.gif#only-dark" alt="Calculating stresses and strains in a phase field crystal model">
    </p>
    <p style="color: var(--md-default-fg-color)">
    Learn how to calculate stresses and strains in the PFC model.
    </p>
</a>
<a href="https://colab.research.google.com/github/vidarsko/ComFiT/blob/main/tutorial/phase_field_crystal_polycrystalline_systems.ipynb" class="card" style="min-width: 160px; flex: 0 1 calc(20.00% - 10px); margin: 5px;">
    <div> <strong> Polycrystalline systems </strong></div>
    <hr>
    <p>
    <img src="img/index_tutorial_pfc_polycrystals.gif#only-light" alt="Polycrystalline phase field crystal simulation">
    <img src="img/index_tutorial_pfc_polycrystals-colorinverted.gif#only-dark" alt="Polycrystalline phase field crystal simulation">
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

## Virtual environment

Using a virtual environment when using ComFiT is highly encouraged because even though we try to write robust code, it is still a library under development, so previously written simulations may break. By keeping your simulations together with the specific version of ComFiT, you make sure that your simulations will not break due to coming updates.

To create a virtual environment, run the following command in your terminal after having navigated to the root folder of your exploration project

```bash
Python -m venv myvenv
```

This will create the folder `myvenv` which will contain the local installation of Python and associated packages.
To activate the virtual environment, simply run

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
