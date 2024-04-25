---
title: 'ComFiT: a Python library for computational field theory with topological defects'
tags:
  - computational field theory
  - topological defects
  - Bose-Einstein condensates
  - phase-field crystal
authors:
 - name: Vidar Skogvoll
   orcid: 0000-0003-4941-6886
   affiliation: 1
 - name: Jonas Rønning
   orcid: 0000-0001-5289-7276
   affiliation: 2
affiliations:
 - name: Department of Physics, University of Oslo, P. O. Box 1048, 0316 Oslo, Norway.
   index: 1
 - name: Okinawa Institute of Science and Technology
   index: 2 
date: 06 July 2023
bibliography: paper.bib
---

# Statement of need

Field theories described by partial differential equations (PDEs) are the backbone of many areas in physics, but simulating them often requires specialized software and programming expertise.
Even though there exist specialized software packages like FEniCS [@alnaesFEniCSProjectVersion2015] for solving PDEs, these are often focused on numerical efficiency at the cost of legibility or user-friendliness and fall short in offering functionalities for visualizing and analyzing outcomes.
In particular, in the realm of many field theories, the study of topological defects - small structures like vortices in fluids - is essential for understanding phenomena such as phase transitions, turbulence and pattern formation.
Due to the shared mathematical structures of these topological defects, recent research has shown that a common computational framework can be used to study them across different physical systems, ranging from Bose-Einstein condensates to nematic liquid crystals and crystalline solids [@skogvollUnifiedFieldTheory2023].
However, a unified computational framework that brings all these systems together is lacking.
ComFiT aims to close this gap, catering to both researchers and educators in physics, by providing a user-friendly, object-oriented framework for setting up a physical system, solving PDEs in one, two and three dimensions, visualizing results, and tracking topological defects.
In so doing, ComFiT also brings advanced models of phase-field crystal modeling and nematic liquid crystals to the Python ecosystem, which are currently scarcely available in the open-source community, especially for three dimensions.

ComFiT sets itself apart from existing open-source Python software for solving PDEs in these ways:

- **Spectral Methods**: Like Dedalus [@burnsDedalusFlexibleFramework2020], ComFiT employs spectral methods for differentiation and integration,  which is more accessible for those familiar with Fourier analysis, unlike more complex finite element/volume approaches common in performance-driven libraries, e.g., FEniCS [@alnaesFEniCSProjectVersion2015], PyClaw [@ketchesonPyClawAccessibleExtensible2012], and Firedrake [@hamFiredrakeUserManual2023].
- **Built-in Visualization**: Unlike the aforementioned libraries, ComFiT includes tailored plotting tools for physical systems, saving users the effort of integrating with external libraries.
- **Topological Defect Analysis**: ComFiT's defect identification and tracking algorithms provide unique insights crucial for studying physical phenomena. To the best of our knowledge, no other integrated Python library for PDEs offers this functionality.

# Summary

The core functionality of ComFiT is provided in the `BaseSystem` class, which defines a computational domain, auxiliary quantities, methods for time evolution, visualization tools and algorithms for identifying and tracking topological defects.
To simulate a physical system, a user can create a class inheriting from `BaseSystem` and implement the specific equations and parameters of the model.
The library also provides a range of predefined models, which inherit from the `BaseSystem` class, such as `QuantumMechanics`, `BoseEinsteinCondensate`, `NematicLiquidCrystal`, and `PhaseFieldCrystal`, each tailored to a specific field theory, which may be used to quickly get started with simulations for research or educational purposes.

![Four example setups of the ComFiT library. (a) The function $f(x,y) = x/a_0-y/a_0$ where $a_0$ is a length scale, (b) a quantum mechanical wavepacket with a nonzero velocity in three dimensions, (c) a Bose-Einstein condensate vortex ring in three dimensions with vortex nodes identified and (d) a square phase-field crystal simulation containing a dislocation dipole. More details of the systems are given in the package documentation, and the code used to make these figures is given in the appendix.](illustration.png)

The project aims to be highly community-driven, continuously improving the code stability, efficiency and usability, and to be a platform for sharing and developing new models and methods for field theories and topological defects.
[The documentation](https://comfitlib.com) is hosted with MkDocs and contains a range of theoretical backgrounds for all the subclasses, tutorials, and examples to guide users through the process of setting up and running simulations, as well as visualizing and analyzing the results.

# Research projects 

ComFiT is a synthesized code base that was originally written in Matlab and used to fuel research into projects exploring the statistical properties and similarities between the dynamics of topological defects between Bose-Einstein condensates and phase-field crystals [@skaugenUnifiedPerspectiveTwodimensional2018] building on the framework for evolving stiff numerical systems described in Ref. [@coxExponentialTimeDifferencing2002] and the methods for tracking defects developed in Refs. [@mazenkoVortexVelocitiesSymmetric1997;@mazenkoVelocityDistributionStrings1999;@anghelutaAnisotropicVelocityStatistics2012].
The code base was further improved, stabilized and extended in the PhD projects of the paper authors, Refs. [@ronningTopologicalDefectsFlows2023;@skogvollSymmetryTopologyCrystal2023].
ComFiT has been used for research into a range of physical systems, including active matter [@ronningSpontaneousFlowsDynamics2023], stirred Bose-Einstein condensates [@ronningClassicalAnalogiesForce2020;@ronningPrecursoryPatternsVortex2023], phase-field crystals [@skogvollDislocationNucleationPhasefield2021;@skogvollStressOrderedSystems2021;@skogvollPhaseFieldCrystal2022;@skogvollHydrodynamicPhaseField2022].
In all of these applications, the role of topological defects is central to understanding the critical dynamics and the statistical properties of the systems, e.g., in Bose-Einstein condensates, the formation of vortices mark the transition to quantum turbulence and in crystals, the formation of dislocations are central to understanding how materials buckle.
The library is at the time of publication in use by the authors in their research on turbulence in 3D active nematics, and dynamics of 3D vortex structures in Bose-Einstein condensates.
The PFC framework is currently under exploration at the Njord Centre at the University of Oslo for modeling crack propagation in Earth's upper crust and the basic framework is used as a foundation for several projects in the theoretical physics group at Porelab, UiO.

# Appendix

The code used to produce the illustration is given below.

```python
!pip install comfit==1.4.2
import comfit as cf
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(10, 10))

# Base System class instance
bs = cf.BaseSystem(2, xlim=[-3,3], ylim=[-3,3])
field = bs.x - bs.y
ax1 = fig.add_subplot(2, 2, 1) 
bs.plot_field(field, ax=ax1)

# # Quantum Mechanical System 
qm = cf.QuantumMechanics(3, xRes=41, yRes=41, zRes=41)
qm.conf_initial_condition_Gaussian(initial_velocity=[0,0.1,0.3])
qm.evolve_schrodinger(200)
ax2 = fig.add_subplot(2, 2, 2, projection='3d')
qm.plot_complex_field(qm.psi, ax=ax2)

# Bose Einstein Condensate System
bec = cf.BoseEinsteinCondensate(3, xRes=41, yRes=41, zRes=41)
bec.conf_initial_condition_Thomas_Fermi()
bec.conf_insert_vortex_ring()
bec.evolve_relax(100)
vortex_nodes = bec.calc_vortex_nodes()
ax3 = fig.add_subplot(2, 2, 3, projection='3d')
bec.plot_field(abs(bec.psi), alpha = 0.2, ax=ax3, colorbar=False)
bec.plot_vortex_nodes(vortex_nodes,ax=ax3)

# Phase-field crystal system 
pfc = cf.PhaseFieldCrystal2DSquare(15,15)
eta = pfc.calc_amplitudes_with_dislocation_dipole()
pfc.conf_PFC_from_amplitudes(eta)
pfc.evolve_PFC(100)
dislocation_nodes = pfc.calc_dislocation_nodes()
ax4 = fig.add_subplot(2, 2, 4)
pfc.plot_field(pfc.psi,ax=ax4)
pfc.plot_dislocation_nodes(dislocation_nodes,ax=ax4,grid=False)
```

# References
