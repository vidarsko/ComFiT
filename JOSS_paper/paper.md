---
title: 'ComFiT: a python library for Computational Field Theory'
tags:
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
Even though there exist specialized software packages for solving PDEs, these are often focused on numerical efficiency at the cost of legibility or user-friendliness and fall short in offering functionalities for visualizing and analyzing outcomes.
In particular, in the realm of many field theories, the study of topological defects - small structures like vortices in fluids - is essential for understanding phenomena such as phase transitions, turbulence and pattern formation.
Due to the shared mathematical structures of these topological defects, recent research has shown that a common computational framework can be used to study them across different physical systems, ranging from Bose-Einstein condensates to nematic liquid crystals and crystalline solids [@skogvollUnifiedFieldTheory2023].
However, a unified computational framework that brings all these systems together is lacking.
ComFiT aims to close this gap, catering to both researchers and educators in physics, by providing a user-friendly, object-oriented framework for setting up a physical system, solving PDEs in one, two and three dimensions, including tools for visualizing results and tracking topological defects.

# Summary

The core functionality of ComFiT is provided in the `BaseSystem` class, which defines a computational domain, auxiliary quantities, methods for time evolution, visualization tools and algorithms for identifying and tracking topological defects.
To simulate a physical system, a user can create a class inheriting from `BaseSystem` and implement the specific equations and parameters of the model.

The library also provides a range of predefined models, which inherit from the `BaseSystem` class, such as `QuantumMechanics`, `BoseEinsteinCondensate`, `NematicLiquidCrystal`, and `PhaseFieldCrystal`, each tailored to a specific field theory, which may be used to quickly get started with simulations or for educational purposes.

The project aims to be highly community-driven, continuously improving the code stability, efficiency and usability, and to be a platform for sharing and developing new models and methods for field theories and topological defects.
The documentation is hosted with MkDocs and contains a range of theoretical backgrounds for all the subclasses, tutorials, and examples to guide users through the process of setting up and running simulations, as well as visualizing and analyzing the results.

# Research projects 
<!-- (Mention (if applicable) a representative set of past or ongoing research projects using the software and recent scholarly publications enabled by it.) -->

ComFiT is a synthesized code base which was originally written in Matlab and used to fuel research into projects exploring the statistical properties and similarities the dynamics of topological defects between Bose-Einstein condensates and phase-field crystals [@skaugenVortexClusteringUniversal2016;@skaugenVelocityStatisticsNonuniform2016;@skaugenOriginInverseEnergy2017;@skaugenSeparationElasticPlastic2018;@skaugenUnifiedPerspectiveTwodimensional2018] building on the framework for evolving stiff numerical systems described in Ref. [@coxExponentialTimeDifferencing2002] and the methods for tracking defects developed in Refs. [@mazenkoVortexVelocitiesSymmetric1997;@mazenkoVelocityDistributionStrings1999;@anghelutaAnisotropicVelocityStatistics2012].
The code base was further improved, stabilized and extended in in the PhD projects of the paper authors, Refs. [@ronningTopologicalDefectsFlows2023;@skogvollSymmetryTopologyCrystal2023].
The framework has fueled research into a range of physical systems, including active matter [@ronningSpontaneousFlowsDynamics2023], stirred Bose-Einstein condensates [@ronningClassicalAnalogiesForce2020;@ronningPrecursoryPatternsVortex2023], phase-field crystals [@skogvollDislocationNucleationPhasefield2021;@skogvollStressOrderedSystems2021;@skogvollPhaseFieldCrystal2022;@skogvollHydrodynamicPhaseField2022].
In all of these applications, the role of topological defects is central to understanding the critical dynamics and the statistical properties of the systems, e.g., in Bose-Einstein condensates, the formation of vortices mark the transition to quantum turbulence and in crystals, the formation of dislocations are central to understanding how materials buckle.
The library is at the time of publication in use by the authors in their research on turbulence in 3D active nematics, and dynamics of 3D vortex structures in Bose-Einstein condensates.
The PFC framework is currently under exploration at the Njord Centre at the University of Oslo for modeling crack propagation in geological upper crust and the basic framework is used as a foundation for several projects in the theoretical physics group at Porelab, UiO.

# References
