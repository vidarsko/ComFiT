# Topological defects

Topological defects are defects that appear due to imposed constraints by boundary conditions[^merminTopologicalTheoryDefects1979].
In this work, we will track topological defects as entities derivable from a coarse defect density field[^skogvollUnifiedFieldTheory2023].
In short, for a system containing topological defects, it is possible to derive a defect density field $\rho$, which upon suitable spatial integration gives the charge of the defects contained in that region.

## The algorithm for tracking topological defects

The algorithm for identifying topological defects is implemented in the method `BaseSystem.calc_defect_nodes`.
It requires an input of a positive real scalar field `defect_density`, which upon integration over a region gives a number proportional to the number of defects in that region.

!!! note "Note"
    The `defect_density` field is not necessarily the same as the defect density field $\rho$ alluded to above.
    Whereas $\rho$ might be a vector- or tensor-valued field, `defect_density` is a real and positive scalar field.
    They are, however, related, and the exact connection depends on the model in question. 
    For example, in the case of a 2D Bose-Einstein condensate, the `defect_density` is given by $|\rho|$, and in three dimension, where $\vec \rho$ is a vector field, `defect_density` is given by $\vec \rho$. 

The algorithm takes`charge_tolerance` and `integration_radius` as inputs, and then follows these steps:

| Step | Illustration|
|----|----|
| 1. Identify the `position_index` corresponding to `max(defect_density)`. | ![Illustration of step 1](img/topological_defects_algorithm_1.png#only-light) ![Illustration of step 1](img/topological_defects_algorithm_1-colorinverted.png#only-dark)  |
| 2. Calculate the integral `charge` of `defect_density` in a ball `region_to_integrate` of radius `integration_radius` around the point corresponding to `position_index`. If `charge>charge_tolerance`, then the point will be added to the identified defect nodes. | ![Illustration of step 2](img/topological_defects_algorithm_2.png#only-light) ![Illustration of step 2](img/topological_defects_algorithm_2-colorinverted.png#only-dark)|
| `while` `charge>charge_tolerance` | |
| &nbsp;&nbsp; 3.1 Store `position_index` in the dictionary `defect_node`.  | |
| &nbsp;&nbsp; 3.2 Calculate position $\mathbf r_0$ of `defect_note` as the expectation value of $(x,y,z)$ by using `defect_density` as a probability distribution function in `region_to_integrate`. | ![Illustration of step 3.2](img/topological_defects_algorithm_3_2.png#only-light) ![Illustration of step 3.2](img/topological_defects_algorithm_3_2-colorinverted.png#only-dark)  |
| &nbsp;&nbsp; 3.3 Add `defect_node` to the list `defect_nodes`. | |
| &nbsp;&nbsp; 3.4 Remove a ball of radius `2*integration_radius` around `position_index` from the region `region_to_search` in which to search for new defect nodes. | ![Illustration of step 3.4](img/topological_defects_algorithm_3_4.png#only-light) ![Illustration of step 3.4](img/topological_defects_algorithm_3_4-colorinverted.png#only-dark) |
| &nbsp;&nbsp; 3.5 Identify the `position_index` corresponding to `max(defect_density)` in the region `region_to_search`. | ![Illustration of step 3.5](img/topological_defects_algorithm_3_5.png#only-light) ![Illustration of step 3.5](img/topological_defects_algorithm_3_5-colorinverted.png#only-dark) |
| &nbsp;&nbsp; 3.6 Calculate the integral `charge` of `defect_density` in a ball `region_to_integrate` of radius `integration_radius` around the point corresponding to `position_index`. | ![Illustration of step 3.6](img/topological_defects_algorithm_3_6.png#only-light) ![Illustration of step 3.6](img/topological_defects_algorithm_3_6-colorinverted.png#only-dark)|

The value of `charge_tolerance` and `integration_radius` must be set depending on the particular system and the nature of defect density.
For example, in the case of a 2D phase-field crystal, the defect density will be $\sqrt{\alpha_{ij} \alpha_{ij}}$, which is the density of Burgers vectors.
Upon integration over a defect node, it is expected to give the absolute value of the Burgers vector of the defect, which is expressed in units of $a_0$, the lattice constant.
Therefore, the input charge tolerance is given in units of $a_0$, e.g., `charge_tolerance=0.2*self.a0`.
In the case of the 2D Bose-Einstein condensate, however, the integral over the defect density will be a unitless number equal to $1$ if integrated over a defect node.
Thus, the input charge tolerance is a unitless number, e.g., `charge_tolerance=0.2`.

Even though it was exemplified in two dimensions, the same algorithm is readily usable in three dimensions, but care needs to be taken to set the proper `charge_tolerance`.
For example, in a 3D Bose-Einstein condensate, the charge density provided is a 2D charge density in 3D space (for details, see Ref.[^skogvollUnifiedFieldTheory2023]).
The integral over the defect density will then no longer be a unitless number, but scale with `integration_radius`, so one might set both `charge_tolerance=0.2*self.a0` and `integration_radius=self.a0` to scale with `a0`.

### Properties of the defect nodes

Typically, `BaseSystem.calc_defect_nodes` is used to make the list of defect nodes, which is then used to calculate the node properties in the model-specific classes.
For instance, `BoseEinsteinCondensate.calc_vortex_nodes` first calls `BaseSystem.calc_defect_nodes` to get the list of defect nodes, and then calculates the properties of the defect nodes, such as their charge.
To calculate the velocity of defects using the method outlined in Ref.[^skogvollUnifiedFieldTheory2023], one must formulate the order parameter as an $n$-component vector field $\vec \psi$ and provide both `psi` (the vector field) and `dt_psi` as an input to `BaseSystem.calc_defect_velocity_field`.
This will return a vector field describing the velocity of the full defect density field, which can be evaluated at the defect nodes to get the velocity of the nodes.
See an example of how this is done in `BoseEinsteinCondensate.calc_vortex_nodes`.




[^merminTopologicalTheoryDefects1979]: Mermin, N. D. (1979). The topological theory of defects in ordered media. Reviews of Modern Physics, 51(3), 591–648. [https://doi.org/10.1103/RevModPhys.51.591](https://doi.org/10.1103/RevModPhys.51.591)

[^skogvollUnifiedFieldTheory2023]: Skogvoll, V., Rønning, J., Salvalaglio, M., & Angheluta, L. (2023). A unified field theory of topological defects and non-linear local excitations. Npj Computational Materials, 9(1), Article 1. [https://doi.org/10.1038/s41524-023-01077-6](https://doi.org/10.1038/s41524-023-01077-6)
