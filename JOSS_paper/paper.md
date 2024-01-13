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
   orcid:
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

Field theories are pivotal in understanding numerous physical phenomena, from quantum mechanics to fluid dynamics. 
Despite their importance, there is a scarcity of comprehensive, user-friendly computational tools catering to this domain. 
Adding to this, in many field theories, topological defects play a crucial role in understanding critical behavior such as phase transitions, but good numerical toolboxes to identify and track these defects are lacking. 
ComFiT addresses this gap, targeting researchers and educators in physics and related fields. 
Its modular design and extensive documentation make it an accessible yet powerful tool for simulating a wide range of field theories.

# Summary

ComFiT adopts an object-oriented approach, centralizing around a BaseSystem class, which provides a foundation for specific models like NematicLiquidCrystal, PhaseFieldCrystal, QuantumMechanics, and BoseEinsteinCondensate. 
Each subclass implements model-specific equations and parameters. 
The library employs Exponential Time Differencing (ETD) methods, particularly ETD2RK and ETD4RK schemes, for accurate time evolution in simulations.
ComFiT has broad applications in academic research, particularly in the study of complex systems in physics. Its capability to simulate diverse field theories makes it a valuable tool for both theoretical and applied research.

# Installation and Documentation
ComFiT can be installed via standard Python package management systems. The library is accompanied by comprehensive documentation, including installation guides, user manuals, and API references, facilitating easy adoption by new users.

# Community and Contributions

ComFiT encourages community involvement. 
The project's repository provides guidelines for code contributions, feature suggestions, and issue reporting, promoting collaborative development and continuous improvement.

# References
