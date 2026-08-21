# Simulation 1: larger and misspecified architectures

This simulation extends the methods-paper study without repeating its matched
architecture design. Data are generated from a four-hidden-node network with
eight Gaussian predictors with an AR(1) correlation of 0.3. Predictors
`x1`--`x5` are active and `x6`--`x8` are exact null covariates. Networks with
2, 4, 8, and 12 hidden nodes
are fitted to Gaussian and Bernoulli responses at sample sizes 500, 1,000, and
2,000.

The primary quantities are meaningful even when the fitted architecture is
wrong:

- grouped Wald rejection for active `x1` and null `x8`;
- convergence and covariance availability;
- delta-method coverage for the response-scale PCE of changing `x1` from 0 to 1;
- prediction RMSE against the true conditional mean.

Individual hidden-layer weights are deliberately not compared across fitted
architectures because hidden-node permutation and architecture mismatch make
such comparisons ill-defined.

The raw generating-network signal is centred and scaled to variance one under
the predictor distribution. The Gaussian response adds independent standard
normal noise, giving a population signal-to-noise ratio of one. The Bernoulli
response uses the same standardised signal as its logit. `x8` is used for the
primary null test because it is the null predictor least correlated with the
active block under the AR(1) design.

## Before running

From the package root, install the exact source being studied and run the tests:

```powershell
R CMD INSTALL .
Rscript -e "devtools::test(reporter = 'summary')"
```

The simulation uses only `nnet`, `statnnet`, and base R. The default full design
contains 24 scenarios, 1,000 replications per scenario, and 10 initialisations
per fitted network. The Hessian is calculated only for the selected fit, avoiding
unnecessary Hessian calculations for discarded initialisations.

Validate the deterministic seeds, exact-null construction, Gaussian and binary
paths, and a repeatability check before starting a long run:

```powershell
Rscript paper/simulations/simulation-1-larger-networks/validate-setup.R
```

## Smoke test

Run two replications for the smallest Gaussian scenario:

```powershell
Rscript paper/simulations/simulation-1-larger-networks/run-simulation.R `
  --run-id=smoke-test --replicates=2 --starts=2 --scenario-id=1 `
  --truth-n=10000 --evaluation-n=200
```

Summarise it with:

```powershell
Rscript paper/simulations/simulation-1-larger-networks/summarise-results.R `
  --run-id=smoke-test
```

## Full run

```powershell
Rscript paper/simulations/simulation-1-larger-networks/run-simulation.R `
  --run-id=main --replicates=1000 --starts=10 --scenario-id=all
```

The runner prints progress and an estimated completion time every ten
replications. It checkpoints every ten replications under `results/main/raw/`.
If interrupted, resume without discarding completed results:

```powershell
Rscript paper/simulations/simulation-1-larger-networks/run-simulation.R `
  --run-id=main --replicates=1000 --starts=10 --scenario-id=all --resume=true
```

An existing run is never replaced silently. `--overwrite=true` first moves it
to a timestamped backup directory.

## Running scenarios in parallel

The safest parallelisation is one process per scenario or per disjoint group of
scenario IDs. Start the first process normally. Once it has printed that the
truth samples are ready, start the other processes with `--resume=true`. For
example, four terminals can run `1,5,9,13,17,21`, `2,6,10,14,18,22`,
`3,7,11,15,19,23`, and `4,8,12,16,20,24` with the same `--run-id=main`.
Each process writes different scenario checkpoint files.

After every process finishes, run:

```powershell
Rscript paper/simulations/simulation-1-larger-networks/summarise-results.R `
  --run-id=main
```

Create the planned six-panel diagnostic figure with:

```powershell
Rscript paper/simulations/simulation-1-larger-networks/plot-results.R `
  --run-id=main
```

## Runtime benchmark

The benchmark runs one replication of every scenario with two starts, separates
fitting from analysis time, and scales fitting time to ten starts and 1,000
replications:

```powershell
Rscript paper/simulations/simulation-1-larger-networks/benchmark-runtime.R
```

Treat the estimate as approximate because optimisation time varies across
replications and machines.

On the development machine, a two-start benchmark covering all 24 scenarios
estimated 26.6 hours for the full 10-start, 1,000-replication design in one R
process. A practical allowance is approximately 30--36 hours sequentially,
8--12 hours with four scenario-level processes, or 5--8 hours with eight
processes. The parallel ranges are deliberately wider than ideal linear scaling
to allow for processor and memory contention. Re-run the benchmark locally
before scheduling the final study.

The full design comprises 24,000 selected networks, 240,000 `nnet`
optimisations, and 24,000 final Hessian/covariance calculations.

## Output contract

Each run records its configuration, scenario grid, PCE truth calculation,
session information, command, checkpoint files, combined results, aggregated
summary, and failure summary. Failed optimisation, augmentation, covariance, or
PCE calculations are retained with their reason and are not silently removed
from denominators.
