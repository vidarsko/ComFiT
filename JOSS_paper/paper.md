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
(A Statement of need section that clearly illustrates the research purpose of the software and places it in the context of related work.)

Field theories are pivotal in understanding numerous physical phenomena, from quantum mechanics to fluid dynamics. 
Despite their importance, there is a scarcity of comprehensive, user-friendly computational tools catering to this domain. 
Adding to this, in many field theories, topological defects play a crucial role in understanding critical behavior such as phase transitions, but good numerical toolboxes to identify and track these defects are lacking. 
ComFiT addresses this gap, targeting researchers and educators in physics and related fields. 
Its modular design and extensive documentation make it an accessible yet powerful tool for simulating a wide range of field theories.

# Summary
(A summary describing the high-level functionality and purpose of the software for a diverse, non-specialist audience.)


ComFiT adopts an object-oriented approach, centralizing around a BaseSystem class, which provides a foundation for specific models like NematicLiquidCrystal, PhaseFieldCrystal, QuantumMechanics, and BoseEinsteinCondensate. 
Each subclass implements model-specific equations and parameters. 
The library employs Exponential Time Differencing (ETD) methods, particularly ETD2RK and ETD4RK schemes, for accurate time evolution in simulations.
ComFiT has broad applications in academic research, particularly in the study of complex systems in physics. Its capability to simulate diverse field theories makes it a valuable tool for both theoretical and applied research.


# Research projects 
(Mention (if applicable) a representative set of past or ongoing research projects using the software and recent scholarly publications enabled by it.)
COMFIT is based on code that has been developed during the authors PhD project and has been used in a range of recent scholarly publications. 
This includes simmulations of active matter in Ref. [^Ronning2023Polar],
stirred Bose-Einstein condensates in Ref. [],
Phasecrystalls in Ref. []
and all of the above in Ref. [^skogvoll2023Topological]

ComFiT encourages community involvement. 
The project's repository provides guidelines for code contributions, feature suggestions, and issue reporting, promoting collaborative development and continuous improvement.

# References
[^Ronning2023Polar]: Rønning, J., Renaud, J., Doostmohammadi, A. & Angheluta, L. (2023). Spontaneous flows and dynamics of full-integer topological defects in polar active matter. Soft Matter, 19, 7513-7527. [https://doi.org/10.1039/D3SM00316G](https://doi.org/10.1039/D3SM00316G)
[^skogvoll2023Topological]: Skogvoll, V., Rønning, J., Salvalaglio, M., Angheluta, L. (2023). A unified field theory of topological defects and non-linear local excitations. npj Comput Mater, 9, 122. [https://doi.org/10.1038/s41524-023-01077-6](https://doi.org/10.1038/s41524-023-01077-6)
