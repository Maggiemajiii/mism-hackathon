# Cell migration — baseline trajectory reconstruction (R)

Turn the migration occupancy grids into cell trajectories with the simplest
principled method: **frame-to-frame Hungarian linking**. This is the *baseline*
every fancier tracker (global LAP, probabilistic/optimal-transport, learned
trackers) has to beat.

The code is deliberately small, dependency-light, and spelled out — you can read
it top to bottom and see exactly how a trajectory is built out of the assignment
problem. Nothing is hidden in a library.

## Files

| file | what it is |
|------|------------|
| `hungarian.R`    | the linear-assignment solver, from scratch (Kuhn's O(n³) potentials method) |
| `tracking.R`     | detect cells → cost matrix → link → chain into trajectories → diagnostics → plot |
| `run_baseline.R` | driver: load every migration movie, print the confidence table, plot the reliable one |

## How it works (one paragraph)

A trajectory does not exist in the data — each frame is just an occupancy grid.
`build_tracks()` manufactures trajectories by, at every frame pair, building the
matrix of distances between the current cells and the next frame's detections,
solving the **assignment problem** (`solve_assignment()`) to match them one-to-one
with least total movement, and then *chaining* each cell's matched index through
all 41 frames. The catch is that this is only trustworthy when a cell's step is
small compared with the spacing between cells — quantified by
`rho = mean step / nearest-neighbour spacing`.

## Run

```bash
Rscript run_baseline.R          # from this folder, or anywhere
```

Requirements: R with **reticulate** and a Python that has **numpy** (identical to
the setup used by `notebooks/data_share.Rmd`, which also reads the `.pkl` files
via `py_load_object`). **clue** is optional — it is used only to speed up the
480-cell movie; without it, the from-scratch solver runs everything (a little
slower on high density).

## Expected output

```
condition  dens   cells   mean_step    nn_space      rho
------------------------------------------------------
no_drug    low        5       3.201      14.142    0.226
no_drug    med       48       2.556       2.236    1.143
no_drug    high     480       1.063       1.000    1.063
with_drug  low        5       1.289       8.062    0.160
with_drug  med       48       1.247       2.914    0.428
with_drug  high     480       0.784       1.000    0.784
```

and `results/tracks_no_drug_low.png`, the reconstructed trajectories of the one
movie we can trust.

Read the table like this: only **no_drug/low** and **with_drug/low & med** have
`rho < 0.5` and are reliably trackable. The drug slows the cells, so their movies
stay trackable to higher density. Everywhere `rho ≳ 1`, the nearest detection is
frequently the wrong cell and the individual tracks are guesses — that is a
property of the data (step ≈ spacing), not of the solver.

## Correctness

`solve_assignment()` was checked against SciPy's `linear_sum_assignment` on 3000
random matrices (identical optimal cost every time), and the full pipeline
reproduces the confidence numbers above.

## This is the baseline — what beats it

Per-frame Hungarian is greedy and memoryless, so it fails exactly where
`rho ≳ 1`. Natural next steps, in order of effort:

1. **Gated cost** — forbid links beyond a few √(2·D·Δt); cheap robustness.
2. **Global / multi-frame LAP** (u-track / TrackMate style) — optimise linking
   over a window so an ambiguous frame is resolved by its neighbours.
3. **Soft / probabilistic association** (entropic optimal transport, MHT) — don't
   commit to one matching; recover *ensemble* motion statistics even when
   individual identities are unrecoverable.
4. **Semi-synthetic validation** — simulate diffusive lattice walks with known
   ground truth to measure link accuracy vs density and pick the method honestly.
