# Changelog

All notable changes to this project will be documented in this file.
Full documentation can be found here: [https://comfitlib.com/](https://comfitlib.com/).

## [Unreleased] - 2.0.0

A cleanup pass over naming consistency and documentation ahead of 2.0.0. No functional/numerical
changes are included in this pass — only renames (breaking, since names are part of the public
API) and documentation fixes, with one exception (see Changed (breaking) below). See ISSUES.md for
consistency issues identified but deliberately deferred (ambiguous naming calls, or changes judged
too large/risky for this pass).

### Renamed (breaking)

- `QuantumMechanics.__init__`: `dimension` parameter renamed to `dim`, for consistency with
  `BoseEinsteinCondensate`, `NematicLiquidCrystal`, and `PhaseFieldCrystal`, which all take `dim`.
- `NematicLiquidCrystal.calc_disclination_nodes_nem` renamed to `calc_disclination_nodes`, for
  consistency with `calc_vortex_nodes` (BEC) and `calc_dislocation_nodes` (PFC), neither of which
  carry a model-specific suffix.
- `NematicLiquidCrystal.calc_disclination_density_nematic` renamed to `calc_disclination_density`,
  for the same reason (also removes a second, differently-spelled redundant-suffix convention that
  coexisted with the `_nem` one above).
- `NematicLiquidCrystal.conf_active_channel`: `d` parameter renamed to `interface_width`, for
  consistency with the same concept in `BoseEinsteinCondensate.conf_dissipative_frame`.
- `BoseEinsteinCondensate.conf_vortex_remover`: `Area` parameter renamed to `area` (lowercase, for
  consistency with every other parameter name in the codebase; PEP8 also reserves capitalized names
  for classes).
- `BoseEinsteinCondensate.evolve_comoving_dGPE`: `velx` parameter renamed to `vel_x`.
- `BoseEinsteinCondensate.calc_harmonic_potential`: `R_tf` parameter renamed to
  `thomas_fermi_radius`, for consistency with the descriptive-parameter-name convention (see
  `docs/Conventions.md`) — `R_tf` was an ad hoc abbreviation of a bare mathematical symbol rather
  than a spelled-out name. This is the same quantity as `QuantumMechanics.conf_harmonic_potential`'s
  `trapping_strength` (both build `1/thomas_fermi_radius**2 * (r - r_mid)**2`), just parametrized
  differently — deliberately, since interactions give the Thomas-Fermi radius a physical meaning
  `QuantumMechanics` lacks, so the two parameters were kept distinct rather than unified.
- `BaseSystem.calc_defect_current_density`: `psi_0` parameter renamed to `psi0`, for consistency
  with the identical equilibrium-amplitude parameter on `calc_defect_density`,
  `calc_defect_density_singular`, and `calc_delta_function`.
- `BaseSystem.get_sym`, `get_sym_tl`, and `get_anti_sym` renamed to
  `get_component_from_symmetric_tensor`, `get_component_from_symmetric_traceless_tensor`, and
  `get_component_from_antisymmetric_tensor` respectively, spelling out what the abbreviated names
  left implicit (see `docs/ClassNematicLiquidCrystal.md`, `AGENTS.md`, and ISSUES.md, all updated
  to match). `get_anti_sym`'s `omega` parameter was also renamed to `tensor` along the way, for
  consistency with the identical parameter on the other two.
- `BaseSystem`: internal `_check_if_fourier_and_adjust` renamed to `check_if_fourier_and_adjust`
  (dropped the leading underscore — it was the only "private" name in the codebase; the codebase
  does not otherwise use underscore-privacy).
- `comfit.tool.tool_math_functions`: `levi_civita_symbol` renamed to `tool_levi_civita_symbol`, for
  consistency with the `tool_` prefix convention used by every other function in `comfit/tool/`.
- `comfit.tool.tool_plotly_colorbar` (internal helpers): `format_tick_value` and
  `generate_numbers_between` renamed to `tool_format_tick_value` and `tool_generate_numbers_between`.
- `comfit.tool`: `tool_add_spacing_2D`/`tool_add_spacing_3D`, `tool_matplotlib_define_2D_plot_ax`/
  `_3D_plot_ax`, and `tool_plotly_define_2D_plot_ax`/`_3D_plot_ax` (and their files) renamed to
  lowercase `2d`/`3d`, per PEP8's lowercase module/function naming convention (now stated
  explicitly in `docs/Conventions.md`) and for consistency with `comfit/phase_field_crystal/`'s
  existing lowercase `2d`/`3d` file naming.
- `QuantumMechanics.conf_initial_condition_Gaussian` renamed to `conf_initial_condition_gaussian`,
  `BoseEinsteinCondensate.conf_initial_condition_Thomas_Fermi` renamed to
  `conf_initial_condition_thomas_fermi`, and `BaseSystem.calc_Gaussian_filter_f`/`calc_Gaussian`
  renamed to `calc_gaussian_filter_f`/`calc_gaussian` — function and variable names are always
  lowercase in this codebase, with no exception for embedded proper nouns (now stated explicitly
  in `docs/Conventions.md`), matching the convention used by NumPy/SciPy/scikit-learn.
- `PhaseFieldCrystal.evolve_PFC_mechanical_equilibrium`: `Delta_t` parameter renamed to `delta_t`
  (PEP8 gives no exception for a capitalized local/parameter name). Confirmed first that this
  quantity is not literally `self.dt` under another name — it's a distinct, coarser interval that
  the method divides `self.dt` into (`round(Delta_t/self.dt)`) to get a step count — and recorded
  the resulting `dt`-vs-`delta_t` naming rule in `docs/Conventions.md` before renaming.
  `NematicLiquidCrystal.calc_dt_psi`'s `delta_t` parameter was already correct under this rule and
  needed no change.
- `NematicLiquidCrystal.calc_disclination_velocity_field`: `T`, `Omega_R`, `g`, `omega` parameters
  renamed to `tangent_vector`, `rotation_vector`, `g_matrix`, `disclination_density_magnitude`.
  This function implements the disclination velocity law of Schimming & Viñals, "Kinematics and
  dynamics of disclination lines in three-dimensional nematics," Proc. R. Soc. A 479:20230042
  (2023) (arXiv:2212.10620), Eq. (8) — confirmed against the paper before renaming, per AGENTS.md's
  requirement that a physics-informed rename be checked rather than guessed. The new names match
  what the call site and the sibling `calc_disclination_density_decoupled`'s docstring already used
  informally. Also filled in the previously-missing `Parameters` docstring entries for these four
  arguments, and fixed `calc_g_matrix`'s `Returns` docstring (said "The disclination density",
  copy-pasted from a neighboring function; now says "The g matrix"). Recorded the general
  descriptive-parameter-naming rule this implies in `docs/Conventions.md`.

### Changed (breaking)

- `BaseSystem.get_component_from_antisymmetric_tensor` (see rename above) took its `tensor`
  argument inconsistently depending on `self.dim`: a length-3 array of component fields in 3D, but
  a bare scalar field (no array wrapper at all) in 2D, flagged by a pre-existing `# TODO` on the
  function (see ISSUES.md). Since a 2D antisymmetric tensor has exactly one independent component,
  it now takes a length-1 array in 2D too, consistent with the 3D convention. Updated the one call
  site, `NematicLiquidCrystal.calc_passive_stress_f` (both the 2D and 3D branches build
  `Antisym_QH`), to match. The computed stress values are unchanged — this only changes the storage
  convention of the intermediate argument.

### Fixed

- `NematicLiquidCrystal.calc_active_force_f` and 17 other call sites across
  `nematic_liquid_crystal.py` (plus one in `BoseEinsteinCondensate.calc_kinetic_energy`) passed a
  generator expression to `np.sum(...)`. Recent NumPy raises `TypeError: Calling np.sum(generator)
  is deprecated` instead of the old (undocumented) behavior of silently falling back to Python's
  builtin `sum()` — which is what these calls actually relied on, since `np.sum(generator)` never
  summed over all axes as the array form would. Replaced with the builtin `sum(...)` directly,
  preserving the original (elementwise-add-the-yielded-arrays) behavior explicitly. This was
  blocking the entire `tests_nematic_liquid_crystal` suite on current NumPy.
- `setup.py`: `kaleido==0.2.1` pinned against an unpinned `plotly` broke on any environment that
  resolved a current `plotly` (7.0 dropped support for Kaleido below v1.0.0 — see
  [plotly.py changelog](https://github.com/plotly/plotly.py)), causing every `plot_save`
  (`fig.write_image`) call to fail with `RuntimeError: Image export requires the Kaleido package,
  v1.0.0 or greater`. Bumped to `kaleido>=1.0.0`, `plotly>=6.1.1`. Kaleido v1 requires a system
  Chrome install rather than bundling one; `.github/workflows/tests_plot.yml` now fetches a
  compatible build via `plotly.io.get_chrome()` before running the plot test suite.

### Documentation

- Large documentation cleanup pass: docstring style/consistency fixes across the codebase, several
  stale/incorrect docstrings corrected to match actual behavior, typo fixes, and removal of an
  orphaned scratch file (`docs/compiled_document.md`). No behavioral changes.
- Follow-up dedicated sweeps deferred from the pass above (see former ISSUES.md entries): a
  package-wide spell-check across every docstring and `docs/*.md` page (more typos found and
  fixed), and a full markdownlint compliance pass across `docs/*.md` (missing image alt text,
  malformed/inconsistent tables, indented code blocks converted to fenced, a duplicate heading, a
  missing top-level heading, and false-positive footnote-reference warnings resolved). Added
  `.markdownlint.json` disabling `MD013` (line-length — doesn't fit this project's prose/math/table
  style) and `MD033` (inline HTML — needed for the card-grid templates in `docs/Templates.md` and
  `docs/index.md`) project-wide, and fixing `MD046` (code-block style) to `fenced`. No behavioral
  changes.

### Structural

- `comfit/tool/tool_plotly_colorbar.py` bundled three distinct functions
  (`tool_plotly_colorbar`, `tool_format_tick_value`, `tool_generate_numbers_between`) in one file,
  unlike the one-function-per-file pattern used everywhere else in `comfit/tool/`. Split into
  `tool_plotly_colorbar.py`, `tool_format_tick_value.py`, and `tool_generate_numbers_between.py`.
  `tool_format_tick_value` and `tool_generate_numbers_between` are now exported directly from
  `comfit.tool` (previously reachable only via `comfit.tool.tool_plotly_colorbar`);
  `tool_plotly_colorbar`'s own import path and behavior are unchanged. Also dropped an unused
  `tool_colormap` import from `tool_plotly_colorbar.py`.

### Reviewed (no change)

- `comfit/plot/plot_vector_field_in_plane_both_plot_libs.py` was flagged as the only file in
  `comfit/plot/` combining both backends in one file instead of following the
  `plot_X_matplotlib.py` / `plot_X_plotly.py` per-backend file-pair pattern. Confirmed with the
  maintainer this was intentional: unlike the other pairs, most of this function's body (marching
  cubes, interpolation, vector scaling) is backend-agnostic setup shared by both backends, with only
  the final rendering calls differing — splitting it would duplicate that shared code across two
  files rather than clean anything up. Left as-is.

- `comfit/nematic_liquid_crystal/plot_field_velocity_and_director_matplotlib.py` and `_plotly.py`
  were flagged as the only backend-split `plot_*` files living outside `comfit/plot/`. Confirmed
  with the maintainer this is intentional: `comfit/plot/` holds generic plotting primitives backing
  `BaseSystem`'s generic `plot_*` methods, while `plot_field_velocity_and_director` is a composite
  visualization specific to nematic physics (overlaying a scalar field, a velocity streamplot, and
  a headless director quiver plotted for the ±n symmetry) that isn't meaningful for other models.
  It is correctly defined only on `NematicLiquidCrystal`, not `BaseSystem`. Left as-is.

## [1.9.6] - 2025-04-24
- Added `calc_coarse_grain` method to the BaseSystem class.

## [1.9.5] - 2025-04-16
- Bug fix for tool_configure_axis function.

## [1.9.4] - 2025-04-16
- Removed the need to pass both `ax` and `fig` to plot functions.

## [1.9.3] - 2025-04-16
- Changed definition of `xlim`. It is now the limit of the x-axis, not the size `xmax` of the domain. See comfitlib.com for more details.
- Added shadows to 3D complex plots with plotly.
- Changed default plotting setting for fourier space from 2pi/a0 to 1/a0. 
- Bug fixes.

## [1.9.2] - 2025-03-25
- Fixed bug related to trying to plot a constant field.

## [1.9.1] - 2025-03-24
- Bug fix: setting ylim with plotly now works for 2D plots. 
- Added `width` and `height` parameters to plot_save-function.
- Added phase angle plot to `plotly` library.
- Changing the requirements to work with numpy>=2.2.4

## [1.9.0] - 2025-03-18
- Added the option to pass `fourier=True` to the plot functions to plot a field that is in Fourier space.

## [1.8.7] - 2025-03-13
- Changed default colormap to 'angle' for plot_complex_field_in_plane function.

## [1.8.6] - 2025-03-13
- Fixed bug in plot_save-function.

## [1.8.5] - 2025-03-11
- Added 'opacity' to 3D plot field functions.

## [1.8.4] - 2025-03-03
- Changed plotting convention cfi.plot_save(n,fig) -> cfi.plot_save(fig,n)

## [1.8.3] - 2025-02-28
- Changed calculation of integrating factors in base system class remove division by zero errors.

## [1.8.2] - 2025-02-21
- Fixed plotly colorbar placement.

## [1.8.1] - 2025-02-21
- Plotly now returns a fig, ax object, like the matplotlib functions.

## [1.8.0] - 2025-02-04
- Created more seamless transition between plotly and matplotlib. 

## [1.7.0] - 2025-01-08
- Changed default plotting library to plotly
- Changes to PFC class

## [1.6.2] - 2024-09-07
- Bug fixes for plotly vector plots.

## [1.6.1] - 2024-09-05
- Added more plotly features for alternative plotting.

## [1.6.0] - 2024-06-30
- Added plotly feature integration for alternative plotting (under development).
- Added polycrystalline methods for the PFC models.

## [1.5.0] - 2024-05-30
- Finalized version of software after JOSS software review. 

## [1.4.2] - 2024-03-18
- Bug fixes and stability improvements.

## [1.4.1] - 2024-03-15
- Fixed error in plot_complex_field function.

## [1.4.0] - 2024-03-15
- Completed writing all the plot functions. 
- Added stress tensor calculations for all the PFC models.
- Completed much of the documentation.

## [1.3.1] - 2024-02-24
- Added stress divergence calculation method to all PFC models.

## [1.3.0] - 2024-02-24
- All mayavi functionality removed and placed in optional extension (ComFiT-mayavi)

## [1.2.2] - 2024-02-24
- Relaxed vtk requirement

## [1.2.1] - 2024-02-24
- Made a wheel for the distribution

## [1.2.0] - 2024-02-24
- Bug fixes and stability improvements.
- Included the Mayavi package for 3D plotting.

## [1.1.0] - 2023-12-07
### Fixed
- Many bug fixes and stability improvements, in particular for the BoseEinsteinCondensate class.

## [1.0.0] - 2023-11-22
- Initial release of the package.
- Containing most functions, but lacking in stability.