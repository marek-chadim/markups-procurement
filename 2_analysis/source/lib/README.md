# 2_analysis/source/lib/ — shared library modules

Utility modules imported by the consumer scripts in `2_analysis/source/`.
Moved here from the `source/` root on 2026-04-15 so the root contains
only pipeline-entry scripts (things `make.sh` invokes directly) while
reusable machinery lives under `lib/`.

## Python modules

### `acf_estimator.py`
Monolithic `ACFEstimator` class for ACF-style production function
estimation (Cobb-Douglas + translog, Kim-Luo-Su 2019 overidentification
lags, Chamberlain 1987 sieve optimal instruments, analytical GMM sandwich
SEs, CWDL measurement-error extensions).

Consumers (10 scripts in `source/`):
- `paper_results.py` — main-text Spec A translog headline
- `adl_instrument_comparison.py` — ABGRS instrument rows 1–9
- `ags_twostep_identification.py` — Andrews-Gentzkow-Shapiro two-step CS
- `beer_c2_table.py` — Beer (2024) C2 replication
- `bmy_czech_analysis.py` — Bond-Hashemi-Kaplan-Zoch decomposition
- `bmy_rolling_stability.py` — rolling-window stability checks
- `cwdl_robustness.py` — Collard-Wexler De Loecker measurement-error robustness
- `dls_markup_comparison.py` — De Loecker-Syverson 8-method comparison
- `dlw_treatment_eval.py` — DLW 2012 treatment evaluation
- `misspecification_diagnostics.py` — translog residual diagnostics

### `dml_core.py`
Shared DML infrastructure: `CrossFitter`, `orthogonal_score_plr`,
`sensitivity_omvb`, `cluster_boot_firm`, `omega_partial_r2`. Used by the
§6.8 ML robustness layer.

Consumers (6 scripts in `source/`):
- `dml_premium.py` — partially linear headline (§6.8.1)
- `dml_sensitivity.py` — Cinelli-Hazlett bias landscape + RV
- `dml_cate.py` — conditional ATE by NACE/reform era
- `dml_iv.py` — partially linear IV (single reform)
- `dml_iv_multi_reform.py` — three-reform IV Anderson-Rubin CS
- `dml_strong_exclusion.py` — ABGRS partial-R² diagnostic

### `style_markups.py`
Matplotlib style helper: Paul Tol bright palette, STIX Two Text serif,
CM mathtext. Call `apply_markups_style()` after matplotlib import.
Consumers: all scripts that produce `.pdf` figures via matplotlib.

## R modules

### `theme_markups.R`
ggplot2 theme helper: Paul Tol palette constants (`markups_blue`,
`markups_pink`), `theme_set(theme_minimal(base_size=11))`, helper
`ggsave_markups()`. Source with `source("lib/theme_markups.R")`.

Consumers: `fect_estimation.R`, `grf_heterogeneity.R`, `lalonde_catt.R`,
`lalonde_overlap.R`, `sunab_event_study.R`.

### `lalonde_functions_helpers.R`
Lalonde (1986) treatment-effect helper functions (propensity reweighting,
CIA checks, balance tables). Source with
`source("lib/lalonde_functions_helpers.R")`.

Consumers: `lalonde_balance.R`, `lalonde_catt.R`, `lalonde_overlap.R`,
`lalonde_placebo.R`, `lalonde_sens.R`.

## Import convention

Python consumers prepend at the top of the file:
```python
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "lib"))
from acf_estimator import ACFEstimator
```

This works regardless of the script's working directory because
`os.path.dirname(__file__)` anchors to the consumer's location, and
`make.sh` does `cd "${MAKE_SCRIPT_DIR}/source"` before invoking scripts.

R consumers change `source("theme_markups.R")` to
`source("lib/theme_markups.R")`. The relative path resolves because
`run_R` in `make.sh` runs from `2_analysis/source/`.
