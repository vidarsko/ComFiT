# Issues

A running list of things identified during ongoing cleanup work (naming consistency,
documentation streamlining, general modernization ahead of/around ComFiT 2.0.0) that were judged
to need more thought, more testing, or a decision from a maintainer before being acted on. Items
here are *not* bugs in the numerics — see AGENTS.md: "Correctness of numerics is the top
priority." This file tracks lower-stakes cleanup debt so it isn't lost. Resolved items should be
removed from this file (their resolution should already be recorded in CHANGELOG.md) rather than
marked done and left in place.

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
