# .legacy/ — historical record of orphaned Stata do-files

This directory holds Stata do-files that were registered in
`master_code.do` at an earlier point but are no longer registered
and no longer produce outputs referenced by the paper. They are
kept on disk as historical record; none of them are executed by
the active pipeline.

## Moved on Apr 14 tidy pass

- `table_panel_didid.do` — staggered-DiD battery (Callaway-Sant'Anna,
  Sun-Abraham, Borusyak-Jaravel-Spiess, Goodman-Bacon, TWFE event
  study) replaced by the PanelIFE AWZ bias-aware section in §5.6
  and the Apr 14 deCdH removal.
- `table_raval_overid.do` — Raval (2023) overidentification test;
  the `app:raval` appendix subsection is disabled via `\iffalse`
  in the paper, so `raval_overid.tex` has no live `\input{}`.
- `table_acf_specification.do` — ACF polynomial/Markov sensitivity;
  `acf_specification.tex` not referenced in the paper.
- `table_aggregate_trends.do` — sales-weighted aggregate markup
  decomposition (Stata port); superseded by the new Stata
  `figure_bmy_decomposition.do` which produces `bmy_decomposition_czech.pdf`
  used in §5.4 + `bmy_decomposition_summary.tex`.

## Why kept on disk

A future referee might request one of these specific robustness
tables, or a replication effort might want to compare. Moving is
fully reversible: `mv .legacy/<file>.do ../` brings any file back
into the active directory, and re-registering in `master_code.do`
takes one line.

See `/Users/marek/.claude/projects/-Users-marek-Desktop-io/memory/matej_mechanism_pivot.md`
for the full restructuring history.
