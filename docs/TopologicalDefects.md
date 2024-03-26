# Topological defects

Topological defects are defects that appear due to imposed constraints by boundary conditions[^merminTopologicalTheoryDefects1979].

In this work, we will track topological defects as course entities derivable from the primary field of observation[^skogvollUnifiedFieldTheory2023].

## Topological defects in a Bose Einstein Condensate

In order to get a sense of what topological defects are, we will look at how they emerge in a Bose-Einstein condensate.
As simulated in the class BoseEinsteinCondensate, the BEC is represented by a complex function $\psi$, which has a specific modulus $\psi_0$ in equilibrium.

The equilibrium condition of $|\psi|=\psi_0$ does not, however, determine the phase $\theta \equiv arg(\psi)$ of $\psi$.
Thus the all states with psi0 are equilibrium states.
The possible equilibrium values of $\theta$ define the _ground state manifold_ $\mathcal M$.
Since the phase of a complex function is an angle between $0$ and $2\pi$, the ground state manifold is the circle.
On top of the ground state manifold, one may define the 

## The algorithm for tracking topological defects

The algorithm for tracking topological defects is bases on the coarse-grained description of topological defects presented in Ref.[^skogvollUnifiedFieldTheory2023].
The algorithm is implemented in the method `BaseSystem.calc_defect_nodes`.
It requires an input of a positive real scalar field `defect_density`, which upon integration over a region gives a number proportional to the number of defects in that region.
The algorithm takes

* `charge_tolerance` 
* `integration_radius`

as inputs, and then follows these steps:

| Step | Illustration|
|----|----|
| 1. Identify the `position_index` corresponding to `max(defect_density)`. | ![Illustration of step 1](img/topological_defects_algorithm_1.png#only-light) ![Illustration of step 1](img/topological_defects_algorithm_1-colorinverted.png#only-dark)  |
| 2. Calculate the integral `charge` of `defect_density` in a ball `region_to_integrate` of radius `integration_radius` around the point corresponding to `position_index`. If `charge>charge_tolerance`, then the point will be added to the identified defect nodes. | ![Illustration of step 2](img/topological_defects_algorithm_2.png#only-light) ![Illustration of step 2](img/topological_defects_algorithm_2-colorinverted.png#only-dark)|
| `while` `charge>charge_tolerance` | |
| &nbsp;&nbsp; 3.1 Store `position_index` in the dictionary `defect_node`.  | |
| &nbsp;&nbsp; 3.2 Calculate position $\mathbf r_0$ of `defect_note` as the expectation value of $(x,y,z)$ by using `defect_density` as a probability distribution function in `region_to_integrate`. | ![Illustration of step 3.2](img/topological_defects_algorithm_3_2.png#only-light) [Illustration of step 3.2](img/topological_defects_algorithm_3_2-colorinverted.png#only-dark)  |
| &nbsp;&nbsp; 3.3 Add `defect_node` to the list `defect_nodes`. | |
| &nbsp;&nbsp; 3.4 Remove a ball of radius `2*integration_radius` around `position_index` from the region `region_to_search` in which to search for new defect nodes. | ![Illustration of step 3.4](img/topological_defects_algorithm_3_4.png#only-light) ![Illustration of step 3.4](img/topological_defects_algorithm_3_4-colorinverted.png#only-dark) |
| &nbsp;&nbsp; 3.5 Identify the `position_index` corresponding to `max(defect_density)` in the region `region_to_search`. | ![Illustration of step 3.5](img/topological_defects_algorithm_3_5.png#only-light) ![Illustration of step 3.5](img/topological_defects_algorithm_3_5-colorinverted.png#only-dark) |
| &nbsp;&nbsp; 3.6 Calculate the integral `charge` of `defect_density` in a ball `region_to_integrate` of radius `integration_radius` around the point corresponding to `position_index`. | ![Illustration of step 3.6](img/topological_defects_algorithm_3_6.png#only-light) ![Illustration of step 3.6](img/topological_defects_algorithm_3_6-colorinverted.png#only-dark)|

The value of `charge_tolerance` and `integration_radius` must be set depending on the particular system and nature of defect density. 
For instance,


[^merminTopologicalTheoryDefects1979]: Mermin, N. D. (1979). The topological theory of defects in ordered media. Reviews of Modern Physics, 51(3), 591–648. [https://doi.org/10.1103/RevModPhys.51.591](https://doi.org/10.1103/RevModPhys.51.591)

[^skogvollUnifiedFieldTheory2023]: Skogvoll, V., Rønning, J., Salvalaglio, M., & Angheluta, L. (2023). A unified field theory of topological defects and non-linear local excitations. Npj Computational Materials, 9(1), Article 1. [https://doi.org/10.1038/s41524-023-01077-6](https://doi.org/10.1038/s41524-023-01077-6)
