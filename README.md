# Synthetic Control with Weight Uncertainty

R code for **Synthetic Control with Weight Uncertainty: Robust Identification
and Statistical Inference**, by Taehyeon Koo and Zijian Guo.

The paper defines the weight-robust treatment effect (WRoTE) through
distributionally robust optimization and develops perturbation-based inference
for it. The repository and main function retain the name `DRoSC`. This code
implements the methods described in the accompanying paper.

## Files

- `src/helpers.R`: `DRoSC()` and supporting functions.
- `src/DRoSC-simulation.R`: simulation-study code.
- `src/Basque.R`: Basque Country analysis.
- `src/semi-real.R`: semi-real-data analyses.

## Requirements

The core method requires `MASS`, `limSolve`, and `intervals`:

```r
install.packages(c("MASS", "limSolve", "intervals"))
```

The empirical and plotting scripts additionally use `Synth`, `ggplot2`,
`ggcorrplot`, `dplyr`, `tidyr`, and `latex2exp`.

## Usage

Run R from the repository root and source the main functions:

```r
source("src/helpers.R")

fit <- DRoSC(
  Y0, Y1, X0, X1,
  lambda = 0,
  Inference = TRUE,
  alpha = 0.05,
  M = 500,
  seed = 1
)

fit$tauHat   # estimated WRoTE
fit$betaHat  # estimated synthetic weights
fit$CI.tau   # perturbation-based confidence set
```

Here, `Y0` and `Y1` are the treated-unit outcomes before and after treatment,
while `X0` and `X1` are the corresponding control-outcome matrices. Each row
of `CI.tau` is one connected component of the aggregated confidence set.
