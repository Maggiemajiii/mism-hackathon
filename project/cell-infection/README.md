# Cell infection — trajectory reconstruction (R)

Runs the **same Hungarian tracking baseline** as `project/cell-migration/`, but on
the infection movies. The tracking code is reused verbatim (`run_infection.R`
sources `../cell-migration/hungarian.R` and `../cell-migration/tracking.R`) — only
the data and the encoding differ.

## Why track the infection movies at all?

The infection question is about *spread*, not motion — but the cells in these
movies still move, and tracking answers three motion questions that feed the
spread analysis:

- **Are the movies reliably trackable?** (the `rho` column)
- **How much do the cells move?** (`um/min`, `%still`) — this tells you whether
  transmission could be carried by moving cells (contact) or happens between
  essentially fixed neighbours (local/free-virus spread).
- **Does the drug change cell motility?** Compare the two conditions — if motility
  is unchanged, any drug effect is on transmission, not movement.

## Encoding note

Infection frames store state (`0` empty, `1` S, `2` I, `3` R), not a plain 0/1
mask. A cell keeps its identity while its state changes and the count is constant
(161), so we track **any occupied site** (`state = NULL` → `value > 0`) with the
exact same linker. (Tracking a single state, e.g. only infected cells, would break
the constant-count assumption, since infected cells appear and disappear.)

## Run

```bash
Rscript run_infection.R
```

Requirements: R + **reticulate** + Python/numpy, and **clue**
(`install.packages("clue")`) — with 161 cells over 401 frames you want the
C-speed assignment solver.

## Expected output

```
condition   exp  cells   mean_step    nn_space      rho    um/min   %still
------------------------------------------------------------------------
no_drug       0    161       0.201       1.414    0.142     0.402      81%
no_drug       1    161       0.197       1.414    0.139     0.394      81%
no_drug       2    161       0.201       1.414    0.142     0.402      81%
with_drug     0    161       0.202       1.414    0.143     0.404      81%
with_drug     1    161       0.200       1.414    0.141     0.399      81%
with_drug     2    161       0.200       1.414    0.141     0.400      81%
```

plus `results/tracks_infection_no_drug_0.png`.

## What it says

`rho ≈ 0.14` everywhere — far below the 0.5 reliability line — so **every**
infection movie is reliably trackable (unlike migration, where medium/high density
failed). The cells are nearly static: mean step ≈ 0.2 site/frame and ~81 % of
frame-to-frame steps are exactly zero, i.e. cells occasionally hop a single
lattice site. Motility (~0.40 µm/min) is essentially identical across drug and
no-drug and across the three replicates, so **the drug acts on transmission, not
on movement** — which is exactly the assumption the SIR/spread analysis relies on.

## Trajectory + state output

`results/trajectory_no_drug_0.csv` is the tidy trajectory table, one row per
(cell, frame): `cell, frame, t_min, x, y, state, state_label`. Tracking supplies
only the **identity** (`cell`); the `state` (0/1/2/3 = empty/S/I/R) is read
verbatim from the frames at each cell's tracked site — ground truth from the data,
not inferred. Every tracked site is occupied, so `state` is always S/I/R.

A useful side effect: a physically valid cell only ever goes S→I→R. If a cell's
`state` jumps *backward* (e.g. R→S), that frame is almost certainly a tracking
mislink — a free correctness check on the linker (~120/161 cells are clean).

This table is the bridge from the motility baseline to the transmission analysis:
with each cell's S→I→R timing and position you can ask *when* and *near whom* every
infection happens.
