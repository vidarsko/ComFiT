# <img src="img/logo.png" width="25" height="25"> ComFiT documentation

ComFiT ([Github](https://github.com/vidarsko/ComFiT)) is a versatile Python library for simulating field theories, including plotting and animation in an object-oriented manner.
If you use ComFiT in your research, please cite the following paper:


!!! quote ""
    Skogvoll, V., & Rønning, J. (2024). ComFiT: A Python library for computational field theory with topological defects. Journal of Open Source Software, 9(98), 6599. [https://doi.org/10.21105/joss.06599](https://doi.org/10.21105/joss.06599)

You can get custom help from the ComFiT assistant by using the GPT model below.
!!! quote ""
    <div style="color: inherit; fill: inherit;"><div style="display: flex;"><div class="notranslate" spellcheck="true" placeholder="Write something, or press 'space' for AI, ' / ' for commands…" contenteditable="true" data-content-editable-leaf="true" style="max-width: 100%; width: 100%; white-space: pre-wrap; word-break: break-word; caret-color: rgb(55, 53, 47); padding: 3px 2px;"><a href="https://chatgpt.com/g/g-6739fa0d28d08191ac2ad0395cd27a6a-comfit-assistant" class="notion-link-mention-token notion-text-mention-token notion-focusable-token notion-enable-hover" data-token-index="0" contenteditable="false" tabindex="0" target="_blank" style="color:inherit;text-decoration:inherit;cursor:pointer" rel="noopener noreferrer"><img style="width:1.2em;height:1.2em;border-radius:3px;vertical-align:-0.15em;margin-right:0.3em" src="https://cdn.oaistatic.com/assets/apple-touch-icon-mz9nytnj.webp"><span class=""><span style="color:rgba(55, 53, 47, 0.65);margin-right:0.3em">ChatGPT</span><span style="border-bottom:0.05em solid solid rgba(55,53,47,.25);font-weight:500;flex-shrink:0">ChatGPT - ComFiT assistant</span></span></a>​</div><div style="position: relative; left: 0px;"></div></div></div>

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
</div>

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
