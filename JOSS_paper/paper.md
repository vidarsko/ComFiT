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
<!-- (A Statement of need section that clearly illustrates the research purpose of the software and places it in the context of related work.) -->

The role of field theories in physics can hardly be understated, serving as a foundational framework for a broad spectrum of physical phenomena.
Despite the availability of numerous technical computational tools designed to tackle partial differential equations (PDEs), these tools often employ complex languages that pose a barrier to individuals not well-versed in computer programming.
Furthermore, these established software packages typically focus on solving PDEs with a particular approach and fall short in offering functionalities for visualizing and analyzing outcomes.
Additionally, in the realm of many field theories, the study of topological defects is essential for understanding phenomena such as phase transitions, turbulence and pattern formation.
However, there is a notable absence of robust numerical toolkits for the identification and monitoring of these defects.

ComFiT aims to close this gap, catering to both researchers and educators in physics.
It introduces an accessible, object-oriented framework for addressing PDEs, with a special emphasis on physical field theories and topological defects.
Rooted in basic Python packages like numpy, matplotlib, and scipy, ComFiT is crafted to be straightforward for users, regardless of their programming proficiency.
The innovative algorithms for depicting and tracking topological defects, integrated into the `BaseSystem` class, provide a cohesive perspective on these defects across a variety of field theories.

ComFiT encourages community involvement.
The project's repository provides guidelines for code contributions, feature suggestions, and issue reporting, promoting collaborative development and continuous improvement.
Particularly well-developed, documented and tested models may be included into the main branch after careful consideration.

# Summary
<!-- (A summary describing the high-level functionality and purpose of the software for a diverse, non-specialist audience.) -->

The core functionality of ComFiT is provided in the `BaseSystem` class, which defines a computational domain, auxiliary quantities, methods for time evolution, visualization tools and algorithms for identifying and tracking topological defects.
The library builds on the exponential time differencing (ETD) methods, particularly ETD2RK and ETD4RK schemes, which employs second and fourth order Runge-Kutta schemes, respectively, for accurate time evolution in simulations.
A user wanting to simulate a physical system can inherit from this class and implement the specific equations and parameters of the model.
The documentation features several tutorials and examples to guide users through the process of setting up and running simulations, as well as visualizing and analyzing the results, which may be used in reserach or in teaching.

The library also provides a range of predefined models, which inherit from the `BaseSystem` class, such as `QuantumMechanics`, `BoseEinsteinCondensate`, `NematicLiquidCrystal`, and `PhaseFieldCrystal`, each tailored to a specific field theory, which may be used to quickly get started with simulations or for educational purposes.
Each subclass implements model-specific equations and parameters, and is well-documented in the ComFiT documentation.

# Research projects 
<!-- (Mention (if applicable) a representative set of past or ongoing research projects using the software and recent scholarly publications enabled by it.) -->

ComFiT is a synthesized code base which was originally written in matlab and used to fuel research into projects exploring the statistical properties and similarities the dynamics of topological defects between Bose-Einstein condensates and phase-field crystals [@skaugenVortexClusteringUniversal2016;@skaugenVelocityStatisticsNonuniform2016;@skaugenOriginInverseEnergy2017;@skaugenSeparationElasticPlastic2018;@skaugenUnifiedPerspectiveTwodimensional2018] building on the framework for evolving stiff numerical systems described in Ref. [@coxExponentialTimeDifferencing2002] and the methods for tracking defects developed in Refs. [@mazenkoVortexVelocitiesSymmetric1997;@mazenkoVelocityDistributionStrings1999;@anghelutaAnisotropicVelocityStatistics2012].
The code base was further improved, stabilized and extended in in the PhD projects of the paper authors, Ref. [@ronningTopologicalDefectsFlows2023;@skogvollSymmetryTopologyCrystal2023].
The framework has fueled research into a range of physical systems, including active matter [@ronningSpontaneousFlowsDynamics2023], stirred Bose-Einstein condensates [@ronningClassicalAnalogiesForce2020;@ronningPrecursoryPatternsVortex2023], phase-field crystals [@skogvollDislocationNucleationPhasefield2021;@skogvollStressOrderedSystems2021;@skogvollPhaseFieldCrystal2022;@skogvollHydrodynamicPhaseField2022], [@skogvollUnifiedFieldTheory2023].
The library is at the time of publication in use by the authors in their research on turbulence in 3D active nematics, and dynamics of 3D vortex structures in Bose-Einstein condensates.
The PFC framework is currently under exploration at the Njord centre at the University of Oslo for modelling crack propagation in upper crust materials and the basic framework is used as a foundation for several projects in the theoretical physics group at Porelab, UiO.

# References
