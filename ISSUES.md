# Issues

A running list of things identified during ongoing cleanup work (naming consistency,
documentation streamlining, general modernization ahead of/around ComFiT 2.0.0) that were judged
to need more thought, more testing, or a decision from a maintainer before being acted on. Items
here are *not* bugs in the numerics — see AGENTS.md: "Correctness of numerics is the top
priority." This file tracks lower-stakes cleanup debt so it isn't lost. Resolved items should be
removed from this file (their resolution should already be recorded in CHANGELOG.md) rather than
marked done and left in place.

## Documentation infrastructure

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
