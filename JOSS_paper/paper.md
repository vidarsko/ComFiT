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

The pivotal role of field theories in physics is profoundly significant, serving as a foundational framework for a broad spectrum of physical phenomena, from the intricacies of quantum mechanics to the complexities of fluid dynamics.
Despite the availability of numerous technical computational tools designed to tackle partial differential equations (PDEs), these tools often employ complex languages that pose a barrier to individuals not well-versed in computer programming.
Furthermore, these established software packages typically focus on solving PDEs with a particular approach and fall short in offering functionalities for visualizing and analyzing outcomes.
Additionally, in the realm of many field theories, the study of topological defects is essential for understanding phenomena such as phase transitions.
However, there is a notable absence of robust numerical toolkits for the identification and monitoring of these defects.

ComFiT emerges as a solution to this void, catering to both researchers and educators in physics and adjacent disciplines.
It introduces an accessible, object-oriented framework for addressing PDEs, with a special emphasis on physical field theories and the nuances of topological defects.
Rooted in basic Python packages like numpy, matplotlib, and scipy, ComFiT is crafted to be straightforward for users, regardless of their programming proficiency.
At its core, the `BaseSystem` class lays the groundwork for engaging with the library and fostering the development of personalized models.
Predefined models such as `NematicLiquidCrystal`, `PhaseFieldCrystal`, `QuantumMechanics`, and `BoseEinsteinCondensate` serve as valuable resources for initiating diverse research endeavors and educational activities.
The innovative algorithms for depicting and tracking topological defects, integrated into the `BaseSystem` class, provide a cohesive perspective on these defects across a variety of field theories.

# Summary
(A summary describing the high-level functionality and purpose of the software for a diverse, non-specialist audience.)

The core functionality of ComFiT is provided in the `BaseSystem` class, which defines a computational domain, auxiliary quantities, methods for time evolution, visualization tools and algorithms for identifying and tracking topological defects.
The library builds on the Exponential Time Differencing (ETD) methods, particularly ETD2RK and ETD4RK schemes, for accurate time evolution in simulations.
A user wanting to simulate a physical system can inherit from this class and implement the specific equations and parameters of the model.
The library also provides a range of predefined models, such as `NematicLiquidCrystal`, `PhaseFieldCrystal`, `QuantumMechanics`, and `BoseEinsteinCondensate`, each tailored to a specific field theory, which may be used to quickly get started with simulations or for educational purposes.
The documentation features several tutorials and examples to guide users through the process of setting up and running simulations, as well as visualizing and analyzing the results, which may be used in reserach or in teaching.

Each subclass implements model-specific equations and parameters. 
The library employs Exponential Time Differencing (ETD) methods, particularly ETD2RK and ETD4RK schemes, for accurate time evolution in simulations.
ComFiT has broad applications in academic research, particularly in the study of complex systems in physics. Its capability to simulate diverse field theories makes it a valuable tool for both theoretical and applied research.


# Research projects 
(Mention (if applicable) a representative set of past or ongoing research projects using the software and recent scholarly publications enabled by it.)
ComFiT is based on code that has been developed during the authors PhD projects, Ref. [^Ronnin2023DefectsAndFlow], [^Skogvoll2023SymmetryTopology].
It have been used in a range of recent scholarly publications,
including simmulations of active matter in Ref. [^Ronning2023Polar],
stirred Bose-Einstein condensates in Ref. [^Ronning2020Classical], [^Ronning2023Precursory],
phase field crystalls in Ref. []
and all of the above in Ref. [^skogvoll2023Topological]. 

ComFiT encourages community involvement. 
The project's repository provides guidelines for code contributions, feature suggestions, and issue reporting, promoting collaborative development and continuous improvement.

# References
[^Ronnin2023DefectsAndFlow]: Rønning, J. (2023). Topological Defects and Flows in BECs and Active Matter. Doctoral Thesis. Univeisity of Oslo [https://www.duo.uio.no/handle/10852/104678](https://www.duo.uio.no/handle/10852/104678)
[^Skogvoll2023SymmetryTopology]: Skogvoll, V. (2023). Symmetry, topology, and crystal deformations: a phase-field crystal approach. Doctoral Thesis. University of Oslo [https://www.duo.uio.no/handle/10852/102731](https://www.duo.uio.no/handle/10852/102731)
[^Ronning2023Polar]: Rønning, J., Renaud, J., Doostmohammadi, A. & Angheluta, L. (2023). Spontaneous flows and dynamics of full-integer topological defects in polar active matter. Soft Matter, 19, 7513-7527. [https://doi.org/10.1039/D3SM00316G](https://doi.org/10.1039/D3SM00316G)
[^Ronning2020Classical]: Rønning, J., Skaugen, A., Hernández-García, E.,  Lopez, C.  & Angheluta L. (2020). Classical analogies for the force acting on an impurity in a Bose–Einstein condensate. New Journal of Physics 22 (7), 073018. [https://doi.org/10.1088/1367-2630/ab95de](https://doi.org/10.1088/1367-2630/ab95de)
[^Ronning2023Precursory]: Rønning, J. & Angheluta, L. (2023). Precursory patterns to vortex nucleation in stirred Bose-Einstein condensates. Physical Review Research 5 (2), 023108 [https://doi.org/10.1103/PhysRevResearch.5.023108](https://doi.org/10.1103/PhysRevResearch.5.023108)
[^skogvoll2023Topological]: Skogvoll, V., Rønning, J., Salvalaglio, M., Angheluta, L. (2023). A unified field theory of topological defects and non-linear local excitations. npj Comput Mater, 9, 122. [https://doi.org/10.1038/s41524-023-01077-6](https://doi.org/10.1038/s41524-023-01077-6)
