# Issues

A running list of things identified during ongoing cleanup work (naming consistency,
documentation streamlining, general modernization ahead of/around ComFiT 2.0.0) that were judged
to need more thought, more testing, or a decision from a maintainer before being acted on. Items
here are *not* bugs in the numerics — see AGENTS.md: "Correctness of numerics is the top
priority." This file tracks lower-stakes cleanup debt so it isn't lost. Resolved items should be
removed from this file (their resolution should already be recorded in CHANGELOG.md) rather than
marked done and left in place.

## Naming — needs a decision, not just an obvious fix

- **`2D`/`3D` vs `2d`/`3d` casing split.** `comfit/tool/` uses uppercase in both filenames and
  function names for six files: `tool_add_spacing_2D.py`/`tool_add_spacing_3D.py`,
  `tool_matplotlib_define_2D_plot_ax.py`/`_3D_plot_ax.py`,
  `tool_plotly_define_2D_plot_ax.py`/`_3D_plot_ax.py`. `comfit/phase_field_crystal/` uses
  lowercase in its six lattice-subclass filenames (`phase_field_crystal_2d_square.py`, etc.).
  These six `tool_*` functions are also part of the public API (re-exported from
  `comfit/__init__.py`) and used at ~80 call sites across `comfit/plot/`,
  `comfit/nematic_liquid_crystal/`, and `comfit/core/base_system_plot.py`. Renaming is mechanical
  but touches a lot of surface area for a purely cosmetic call; needs a decision on which casing
  wins before it's worth doing as one coordinated change (files + `__init__.py` exports + every
  call site + docs).

- **`Gaussian`/`Thomas_Fermi` vs `ordered`/`disordered` capitalization.** Method names like
  `conf_initial_condition_Gaussian`, `conf_initial_condition_Thomas_Fermi`, `calc_Gaussian_filter_f`,
  `calc_Gaussian` capitalize the embedded distribution name, while `conf_initial_condition_ordered`,
  `conf_initial_condition_disordered` don't (because "ordered" isn't a proper noun). This might
  already be the "right" convention (capitalize eponyms, lowercase ordinary adjectives) rather than
  an inconsistency — flagging for a maintainer call rather than assuming either direction.

- **`calc_disclination_velocity_field`'s parameter names** (`T`, `Omega_R`, `g`, `omega`) are
  terse/cryptic compared to the descriptive-name convention used almost everywhere else in
  `NematicLiquidCrystal` (e.g. `dipole_vector`, `charge_tolerance`). Renaming them to something like
  `tangent_vector`, `rotation_vector`, `g_matrix`, `omega_field` seems right but requires someone
  who understands the physics to confirm those are the correct descriptive names before it's done
  as a "safe" rename — left alone this pass.

- **`dt`/`delta_t`/`Delta_t` for "time step".** `BaseSystem` sets `self.dt` (default `0.1`) in
  `base_system_init.py`. `NematicLiquidCrystal` evolvers use a local `delta_t` in places,
  `PhaseFieldCrystal` evolvers use `Delta_t` (also a PEP8 violation — capitalized local variable) in
  others. Before renaming, need to check whether these locals are literally `self.dt` under another
  name or a distinct per-substep quantity (e.g. an RK stage sub-step) — conflating them would be a
  behavior-relevant mistake, not just a rename.

- **`R_tf` (BoseEinsteinCondensate harmonic trap) vs `trapping_strength` (QuantumMechanics harmonic
  trap).** These parametrize conceptually related harmonic traps differently (radius vs. strength),
  so this is likely not actually the "same quantity" and may not need reconciling at all — noted so
  it isn't silently swept into a batch rename later without checking.

## Structural / organizational (bigger than a rename)

- **`comfit/plot/plot_vector_field_in_plane_both_plot_libs.py`** is the only file in `comfit/plot/`
  that doesn't follow the `plot_X_matplotlib.py` / `plot_X_plotly.py` per-backend file-pair pattern
  used by the other 9 matched pairs in that directory. Worth considering whether to split it to
  match, but that's a structural change, not a rename — deferred.

- **`comfit/nematic_liquid_crystal/plot_field_velocity_and_director_matplotlib.py` and
  `_plotly.py`** are the only backend-split `plot_*` files that live outside `comfit/plot/` (every
  other model dispatches through `comfit/plot/` via `BaseSystem.plot_*`). Possibly intentional
  (model-specific plot, not a generic `BaseSystem` one) — flagging as an architecture question
  rather than assuming they should move.

- **`comfit/tool/tool_plotly_colorbar.py`** now correctly has all-`tool_`-prefixed function names
  (`tool_plotly_colorbar`, `tool_format_tick_value`, `tool_generate_numbers_between`), but still
  bundles three distinct utilities in one file, unlike the one-function-per-file pattern everywhere
  else in `comfit/tool/`. Worth splitting, but a file split changes import paths too — deferred as
  a judgment call on whether it's worth the churn for three small, tightly related functions.

## Documentation infrastructure

- **`docs/compiled_document.md`** looks like an orphaned generated/scratch file — it's tracked in
  git, but not referenced from `mkdocs.yml`'s nav, and its (only) content is a near-duplicate
  fragment of `docs/ClassBaseSystem.md`'s ETD2RK section. Probably safe to delete, but leaving it
  for a maintainer to confirm rather than deleting tracked content unilaterally.

- **Repo-wide markdownlint violations.** `docs/Conventions.md` mandates markdownlint for all
  Markdown docs, but the IDE's markdownlint integration flags many pre-existing violations across
  `docs/*.md` just from files touched in this pass (trailing whitespace, `MD060` table-pipe-spacing,
  `MD033` inline HTML in Templates.md-style card blocks, `MD012` multiple consecutive blank lines).
  This wasn't a targeted lint pass — a dedicated sweep (run markdownlint across all of `docs/`, fix
  what it flags, decide whether the inline-HTML card templates need a `.markdownlint.json` override
  instead of "fixing") is its own follow-up task, not something to fix opportunistically file by
  file.

- **Full package-wide spell-check.** Typos were fixed opportunistically wherever this pass already
  had a file open for another reason (see CHANGELOG.md's "Documentation" section for the ones
  found). A dedicated spell-check pass across every docstring and every `docs/*.md` page has not
  been done and would likely turn up more.

## Design smells noticed in passing (functional, not naming — flagged for awareness only)

- `BaseSystem.get_anti_sym` has a pre-existing `# TODO: I don't like that the input vector is a
  scalar field in 2 dimensions. (Vidar 11.03.24)` comment — a real API-shape concern (not just a
  name), already flagged by the author previously. Left untouched; out of scope for a
  naming/documentation pass, and touching it would be a functional change.
