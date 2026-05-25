# Law Review Submission System: Agent-Based Simulation

Replication code for:

> Topaz, C. M. (2026). Inefficiency and inequity of the law review submission system. *PLOS ONE*.

## Overview

This repository contains the R simulation code for an agent-based model of the law review submission market. The model compares the current decentralized expedite system against centralized matching (deferred acceptance) and several partial reforms, across dimensions of efficiency (total match quality) and equity (credential-based advantage).

## Contents

| File | Purpose |
|---|---|
| `simulation.R` | Core simulation: defines all eight mechanisms, runs the baseline Monte Carlo, the four sensitivity analyses, the 399-cell heatmap grid, and the convergence diagnostics. Saves figures and `.rds` outputs. |
| `regenerate_figures.R` | Regenerates all ten manuscript figures from saved `.rds` files without re-running any simulations. Useful for adjusting figure styling without paying the full Monte Carlo cost. |
| `triage_robustness.R` | Editorial-triage robustness check added during peer review (PLOS ONE revision). Implements an explicit desk-rejection step variant and runs a focused Monte Carlo comparison (Current vs. Current + Triage vs. DA). |

## Requirements

- R (version 4.0 or later recommended; tested with 4.5.1)
- Required packages (installed automatically if missing): `parallel`, `pbapply`, `data.table`, `ggplot2`, `scales`, `patchwork`

## Usage

### Full pipeline (baseline + all sensitivities + heatmap; ~1–2 hours on a multi-core machine)

```r
RUN_PIPELINE <- TRUE
source("simulation.R")
```

### Baseline Monte Carlo only (~10 minutes)

```r
RUN_PIPELINE <- "baseline"
source("simulation.R")
```

### Load function definitions without running anything

```r
source("simulation.R")  # by default RUN_PIPELINE is FALSE
```

### Triage robustness check (~5 minutes)

```r
source("triage_robustness.R")
```

### Regenerate figures from saved `.rds` files

```r
source("regenerate_figures.R")
```

## Output

Numerical results are written to `output/` as `.rds` files; figures are written to `manuscript/plos/figures/` as PLOS-compliant TIFFs (5.2 in width, 300 dpi, LZW compression). Both directories are created automatically.

## Reproducibility

All Monte Carlo runs use fixed seeds (`set.seed(2026)` and `clusterSetRNGStream(cl, iseed = 2026)`) for full reproducibility. Within each replication, all mechanisms operate on the same population of journals, articles, editorial assessments, queue orders, and random draws (common random numbers), so mechanism comparisons reflect mechanism design rather than sampling variation.

## Citation

If you use this code, please cite:

```bibtex
@article{topaz2026lawreview,
  author  = {Topaz, Chad M.},
  title   = {Inefficiency and inequity of the law review submission system},
  journal = {PLOS ONE},
  year    = {2026},
  doi     = {[to be added]}
}
```

## License

MIT License. See [LICENSE](LICENSE) for details.

## Contact

Chad M. Topaz · cmt6@williams.edu
