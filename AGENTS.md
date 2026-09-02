# Agent instructions for ComFiT

ComFiT is a published, PyPI-distributed research library for computational
field theory (Schrödinger equation, damped Gross-Pitaevskii, phase-field
crystal, active nematics, etc.), used by researchers to produce results they
publish on. Correctness of numerics is the top priority — more important than
style, coverage of edge cases, or speed of changes.

## Coding style

Follow `docs/Conventions.md`: PEP8, NumPy-style docstrings (a template for
both stand-alone functions and model methods is given there), a fixed import
order (stdlib, third-party, local — blank line between groups), and fixed
notational choices (e.g. `\mathfrak i` for the imaginary unit, `_f` suffix
for Fourier-space quantities). Match existing docstrings rather than
inventing a new style.

## Documentation writing style

Derived from reading `docs/*.md` (in particular `ClassBaseSystem.md`,
`Plotting.md`, `index.md`, `TopologicalDefects.md`) and validated by
drafting sample sections against them (see git history of this file for the
trial-and-error). Match this when writing or editing docs pages — it's a
distinct register from code docstrings (covered above) or this AGENTS.md
file itself.

**Pick the register by what the content is, not by habit.** Four show up;
don't default to derivation for everything:

1. **Derivation** (`ClassBaseSystem.md`'s ETD2RK/Coarse-graining sections) —
   for a result the reader could plausibly want to re-derive or check: walk
   the algebra step by step, one `$$...$$` block per step, before stating
   the final result.
2. **Procedure** (`TopologicalDefects.md`'s node-finding algorithm) — for an
   algorithm the reader will trace through rather than re-derive: a
   numbered step list as a Markdown table with an illustration image in the
   second column per step. Nested/looped steps get decimal numbering
   (`3.1`, `3.2`, ...) with `&nbsp;&nbsp;` indent inside the same cell
   rather than a nested list.
3. **Reference** (`Plotting.md`'s keyword tables) — for a parameter with one
   fixed, stable meaning (constructor kwargs, plot kwargs): a flat
   `| Keyword | Definition | Default value |` table. When a default depends
   on `dim`, express it as a stacked brace matched to dim 1/2/3, not a
   sentence of "if...else".
4. **Justification-by-worked-case** (`TopologicalDefects.md`'s
   `charge_tolerance`/`integration_radius` discussion) — for reasoning that
   depends on the model: prose, one paragraph per concrete case, ending in
   an inline code example (`charge_tolerance=0.2*self.a0`). Don't reach for
   a table just because the content has a tabular shape (dimension × value)
   — if the point is *why* a value is what it is, tabulating it deletes the
   reasoning that was the actual point.

A standard/definitional result (a determinant, a normalization constant)
gets stated and cited, not force-derived — padding a page with algebra
nobody asked for is as wrong as skipping algebra the reader needs.

**Cross-cutting rules, regardless of register:**

- **Open in medias res, not with a definition.** Real sections start
  mid-claim ("A general feature that will be reused is that of vortex
  fields."), not a textbook opening ("A topological defect is a point
  where..."). Define a term the first time it's *used* in a formula or
  claim, not in a preceding throat-clearing sentence.
- **Ground an abstract formula in per-model concrete instances.** Right
  after stating a general result, give one or two sentences showing what it
  reduces to for specific models ("in the case of a 2D Bose-Einstein
  condensate, the `defect_density` is given by $|\rho|$...").
- **"We" for derivations, "you" for instructions.** Math/theory sections use
  first-person plural ("we get", "we therefore have to evaluate..."). Practical
  how-to sections (installation, contributing) switch to direct imperative
  address ("Fork the repository", "You can install ComFiT using...").
- **Tie every symbol to its code name, immediately.** The moment a
  mathematical object is defined, name the exact attribute/function that
  implements it in backticks (`bs.dif[0]`, `calc_wavenums`,
  `calc_nonlinear_evolution_function_f`). Equations and code are interleaved,
  not siloed into separate "theory" and "usage" sections.
- **Use mkdocs-material admonitions with purpose, not decoration:**
  `!!! equation "Title"` for a boxed key result, `!!! note "Title"` for a
  clarifying aside, `??? example "Title"` / `??? abstract "Title"` for a
  collapsible worked example, `=== "tab"` blocks for parallel
  matplotlib/plotly variants. Give each admonition a specific, descriptive
  title, not a generic one.
- **Figures** get both a `#only-light` and `#only-dark` image variant back to
  back, followed by an italicized caption in the form `*Short label:*
  longer description.`
- **Cite with footnotes**: `[^somekey]` inline, full reference collected at
  the bottom of the file under a matching `[^somekey]:` entry — see the
  bottom of `ClassBaseSystem.md`.
- **Version-pin snippets that quote internal source** (not public API) with
  a leading comment, e.g. `# In ComFiT 1.9.0`, since that code will drift —
  this tells the reader the snippet is a point-in-time illustration, not a
  live guarantee.
- **Close a subsection with a one-line forward pointer** when a subtlety is
  deliberately deferred ("We will cover this next.", "We will detail this
  next.") rather than leaving the reader to wonder if it was forgotten.
- **Plain, direct sentences; don't over-hedge.** It's fine to leave a stub
  section with just a heading and "To be written" rather than either
  skipping it or padding it out.
- Write with correct spelling and grammar even where the tone stays
  informal — don't carry forward incidental typos from existing pages when
  editing them.

## Library overview

Checked against the source (`comfit/core`, `comfit/quantum_mechanics`,
`comfit/bose_einstein_condensate`, `comfit/nematic_liquid_crystal`,
`comfit/phase_field_crystal`, `comfit/plot`, `comfit/tool`) as of this
writing. Method signatures and behavior can still change — grep for a
symbol before relying on it rather than trusting this list blindly, and
update this section when they drift.

ComFiT (`import comfit as cf`) simulates field theories on grids with
periodic boundary conditions.

**Class hierarchy.** `BaseSystem` (`comfit/core/base_system.py`, itself
assembled from `BaseSystemInit/Conf/Evolve/Calc/Plot/Get` mixins) sets up
the grid and provides shared calc/plot machinery but has no dynamics or
`conf_*` methods of its own (`BaseSystemConf` is literally empty — all
`conf_*` methods live on the model subclasses). Four models inherit from it:
`QuantumMechanics`, `BoseEinsteinCondensate`, `NematicLiquidCrystal`,
`PhaseFieldCrystal`. Constructor signatures are **not uniform**:
`BoseEinsteinCondensate(dim, **kwargs)` and `NematicLiquidCrystal(dim,
**kwargs)` take `dim` first; `QuantumMechanics(dimension, **kwargs)` names
the same argument `dimension` instead. `PhaseFieldCrystal` itself is not
meant to be instantiated directly — users instantiate one of six concrete
subclasses exported from `comfit/__init__.py`
(`PhaseFieldCrystal1DPeriodic(nx, **kwargs)`,
`PhaseFieldCrystal2DTriangular(nx, ny, **kwargs)`,
`PhaseFieldCrystal2DSquare(nx, ny, **kwargs)`,
`PhaseFieldCrystal3DBodyCenteredCubic(nx, ny, nz, **kwargs)`,
`PhaseFieldCrystal3DFaceCenteredCubic(nx, ny, nz, **kwargs)`,
`PhaseFieldCrystal3DSimpleCubic(nx, ny, nz, **kwargs)`), which take
unit-cell counts rather than `dim`/`xRes` and derive the grid from the
crystal's lattice.

**Config/state attributes** (set in `BaseSystemInit.__init__`): `dim`
(1/2/3), `dx`/`xmin`/`xmax`/`xlim`/`xRes` (and y/z equivalents when
`dim>1`), `dt` (default `0.1`), `plot_lib` (`'plotly'` or `'matplotlib'`,
**default `'plotly'`**). Derived: `psi`/`psi_f` (primary field and its
Fourier transform — name/shape vary by model, e.g. `Q` for the nematic
tensor), `x`/`xmid`/`xmidi`/`size_x` (and y/z), `Res`, `dims`, `rmin`/`rmax`
(always length-3 lists `[xmin,ymin,zmin]`/`[xmax,ymax,zmax]`, regardless of
`dim`), `volume`, `dV`, `time`, `k` (wavenumbers per axis), `dif`
(`1j*k[i]`, for differentiation), `a0` (length scale, default `1`).
**Gotcha:** `rmid` is *not* consistently shaped like `rmin`/`rmax` — it's
the scalar `xmid` when `dim==1` but `[xmid, ymid]` / `[xmid, ymid, zmid]`
when `dim` is 2/3. Code that handles `rmid` generically across dimensions
needs to special-case `dim==1`.

**Broadcasting:** coordinate arrays are shaped for broadcasting without
`meshgrid` — e.g. in 2D, `x.shape == (xRes, 1)` and `y.shape == (1, yRes)`,
so `x + y` broadcasts to `(xRes, yRes)`.

**Method naming convention:** `calc_*` computes and returns a value;
`conf_*` mutates state (sets `psi`/`psi_f` or similar), returns `None`, and
only exists on model subclasses, not `BaseSystem`; `evolve_*` advances the
simulation state (increments `time` by `dt` per step), returns `None`;
`plot_*` returns `(fig, ax)`; `get_*` — narrower than the name suggests, in
practice only `get_sym`/`get_sym_tl`/`get_anti_sym` on `BaseSystem`, for
(anti)symmetrizing tensor components, not a general accessor pattern.

**Fourier fields** are suffixed `_f`. Transform via `self.fft`/`self.ifft`
(`BaseSystem.fft`/`.ifft`, both in `base_system.py`), which call
`scipy.fft.fftn`/`ifftn` **restricted to the last `self.dim` axes**
(`axes=range(-self.dim, 0)`). This means a field can carry extra leading
axes (e.g. tensor/vector components, as with the nematic `Q` tensor) that
`fft`/`ifft` leave alone and only the trailing spatial axes get
transformed. Derivatives: `self.ifft(self.dif[0] * field_f)` for
∂/∂x, Laplacian via `self.ifft(-self.calc_k2() * field_f)` (`.real` if the
field is real).

**Time evolution:** `BaseSystem` itself has none; each model implements its
own (`evolve_schrodinger`, `evolve_dGPE`/`evolve_relax`/
`evolve_comoving_dGPE`, `evolve_nematic`/`evolve_nematic_no_flow`,
`evolve_PFC`/`evolve_PFC_mechanical_equilibrium`/
`evolve_PFC_hydrodynamic`). A custom model subclasses `BaseSystem`, defines
a linear operator (`omega_f`) and a nonlinear evolution function, and
drives them through `self.calc_integrating_factors_f_and_solver(omega_f,
method='ETD2RK')` (also `'ETD4RK'`, both implemented in
`base_system_calc.py`/`base_system_evolve.py`) inside a hand-written
`evolve` loop.

**Plotting:** `plot_field`, `plot_complex_field`, `plot_angle_field`,
`plot_vector_field`, their `*_in_plane` variants, and `plot_nodes` are
defined **once**, on `BaseSystem` (`base_system_plot.py`), and dispatch
internally to a `*_matplotlib`/`*_plotly` free function in `comfit/plot/`
based on the `plot_lib` kwarg (falls back to `self.plot_lib`) — there is no
per-model override needed for the common cases. `plot_nodes(nodes, ...)` is
generic across defect types: it's the same method used for vortex nodes
(BEC), disclination nodes (nematic), and dislocation nodes (PFC) — the
`nodes` argument is whatever `list[dict]` the model's own
`calc_*_nodes`-style method produced. Models do add their own plot methods
on top: `PhaseFieldCrystal` overrides `plot_field` and adds `plot_PFC`,
`plot_orientation_field`; `NematicLiquidCrystal` adds
`plot_field_velocity_and_director`. `plot_subplots(rows, cols)` returns
`(fig, axs)` — `axs` is a flat list if either dimension is 1, else a 2D
array; individual `plot_*` calls take `fig=`/`ax=` to target a subplot.
Animations: loop evolve + plot + `self.plot_save(fig, n)`, then
`cf.tool_make_animation_gif(number_of_frames - 1)` (or
`tool_make_animation_movie`).

**Per-model notables** (verified method names, not exhaustive):

- `QuantumMechanics`: `evolve_schrodinger`, `conf_initial_condition_Gaussian(position, width, initial_velocity)`,
  `conf_harmonic_potential`, `conf_hydrogen_state`, `conf_wavefunction`,
  `calc_hydrogen_state`.
- `BoseEinsteinCondensate`: `evolve_dGPE`, `evolve_relax`,
  `evolve_comoving_dGPE`, `conf_initial_condition_disordered`,
  `conf_initial_condition_Thomas_Fermi`, `conf_external_potential`,
  `conf_insert_vortex`, `conf_insert_vortex_dipole`,
  `conf_insert_vortex_filament`, `conf_insert_vortex_ring`,
  `conf_vortex_remover`, `conf_dissipative_frame(interface_width=7, ...)`,
  `calc_vortex_nodes` (fed into the generic `plot_nodes`), `calc_velocity`,
  `calc_superfluid_current`, `calc_hamiltonian`.
- `NematicLiquidCrystal`: field is the tensor `Q`, evolved via
  `evolve_nematic`/`evolve_nematic_no_flow`;
  `conf_initial_condition_ordered`, `conf_insert_disclination_dipole`,
  `conf_initial_disclination_lines`, `conf_active_channel`, `conf_velocity`,
  `calc_active_force_f`, `calc_passive_force_f`, `calc_pressure_f`,
  `calc_disclination_density_nematic`, `calc_order_and_director`,
  `calc_disclination_nodes_nem` (note the `_nem` suffix — not
  `calc_disclination_nodes`).
- `PhaseFieldCrystal` (base, shared by all six lattice subclasses): field
  `psi` is a real scalar representing crystalline density, evolved via
  `evolve_PFC`; `conf_PFC_from_amplitudes`/`calc_PFC_from_amplitudes`,
  `conf_apply_distortion`, `conf_strain_to_equilibrium`,
  `conf_create_polycrystal`,
  `calc_nonlinear_evolution_function_conserved_f`/`_unconserved_f`,
  `calc_dislocation_density`, `calc_dislocation_nodes`,
  `calc_orientation_field`, `calc_free_energy`, `calc_stress_tensor`.
