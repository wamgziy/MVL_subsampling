# Simulation Data for "Optimal Subsampling for Multivariate Longitudinal Data"

This repository contains the R script used to generate the simulated data
sets in the simulation study (Section 4) of

> **Optimal Subsampling for Multivariate Longitudinal Data**
> Ziyang Wang, HaiYing Wang, Nalini Ravishanker.

The script produces the four full-data sets that correspond to the four
scenarios in Figure 1 of the paper. Subsampling and estimation code is
not included here.

## What the script produces

`simu.R` writes four `.rds` files to the working directory:

| File             | Covariate dist. | n_i (≈) | n_i probabilities                |
|------------------|-----------------|---------|-----------------------------------|
| `data_n5.rds`    | MVN             | 5       | {4,5,6} with (1/20, 9/10, 1/20)   |
| `data_n10.rds`   | MVN             | 10      | {8,9,10} with (1/20, 1/20, 9/10)  |
| `data_n5_t.rds`  | t with df = 2   | 5       | {4,5,6} with (1/20, 9/10, 1/20)   |
| `data_n10_t.rds` | t with df = 2   | 10      | {8,9,10} with (1/20, 1/20, 9/10)  |

Each file is a list `(data, nii)`:

- `data` is a long-format matrix with one row per (subject, time, attribute).
  The first 13 columns are named: `id`, `int`, `X1`–`X4`, `D1`, `D2`, `Y`,
  `int1`, `time1`, `int2`, `time2`. The remaining (unnamed) columns hold
  the block-diagonal fixed-effects design `X^v_i = I_L ⊗ X_i`.
- `nii` is the integer vector of per-subject observation counts.

## Model

For subject *i* with *n_i* time points and *L* response attributes, the
linear mixed-effects model is

```
Y_i = X_i B + Z_i Γ_i + E_i,           i = 1, ..., m,
```

or in vectorized form

```
y^v_i = X^v_i β^v + Z^v_i γ^v_i + ε^v_i,
```

with `X^v_i = I_L ⊗ X_i`, `Z^v_i = I_L ⊗ Z_i`,
`γ^v_i ~ N_{qL}(0, D)`, and `ε^v_i ~ N_{n_i L}(0, Σ ⊗ I_{n_i})`.
See Section 2 of the paper for full details.

## Simulation settings (Section 4)

- `m = 500,000` subjects, `L = 2` attributes
- Fixed-effects vector `β^v = (5, 10, 15, ..., 50)'` (length 10; with
  `L = 2` this gives `d = p = 5` fixed effects per attribute)
- `q = 2` random effects per attribute (random intercept + random slope),
  so `D` is `(qL) × (qL) = 4 × 4`:

  ```
  D = [[4,   1,   1.5, 2  ],
       [1,   5,   2,   2.5],
       [1.5, 2,   3,   2  ],
       [2,   2.5, 2,   4  ]]
  ```

- Within-subject error covariance: `Σ ⊗ I_{n_i}` with `Σ = diag(3, 2)`
- Covariate scenarios: rows of `X_i` from MVN or multivariate t with
  `df = 2`; the `(p − 1)` non-intercept covariates use a compound-symmetric
  cross-column correlation (φ = 0.5) and are independent across time

## Requirements

- R (≥ 4.0)
- R packages: `MASS`, `mvtnorm`, `doParallel`

Install the packages with

```r
install.packages(c("MASS", "mvtnorm", "doParallel"))
```

## Usage

```bash
Rscript simu.R
```

Each scenario generates `m = 500,000` subjects and `L = 2` attributes per
subject. Generation is parallelized via `doParallel`; by default it uses
`detectCores() − 1` workers. The output files are large — total runtime
and disk usage scale with the number of cores and available memory.
