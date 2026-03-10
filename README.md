# Law Review Submission System: Agent-Based Simulation

Simulation code for the paper:

> Topaz, C. M. (2026). Inefficiency and inequity of the law review submission system. *PLOS ONE*.

## Overview

This repository contains the R simulation code for an agent-based model of the law review submission market. The model compares eight article-placement mechanisms — including the current decentralized expedite system and centralized alternatives such as deferred acceptance — across dimensions of efficiency (total match quality) and equity (credential-based advantage).

## Requirements

- R (version 4.0 or later recommended)
- Required packages (installed automatically if missing): `parallel`, `pbapply`, `data.table`, `ggplot2`, `scales`, `patchwork`

## Usage

```r
source("simulation.R")
```

The script will install any missing packages, run the simulation, and save output figures and data files. By default, figures are saved to `manuscript/figures/` and numerical output to `output/`. These paths can be changed by editing `FIGURE_DIR` and `OUTPUT_DIR` at the top of the script.

## License

This code is made available under the MIT License. See [LICENSE](LICENSE) for details.
