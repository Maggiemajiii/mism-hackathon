# Data

Raw movies (not tracked in git; see `.gitignore`) live in `data/raw/` as pickled
dicts `{"A_list": ndarray (T, X, Y), "metadata": {...}}`.

| group | files | frames | interval | cells/field | encoding |
|-------|-------|--------|----------|-------------|----------|
| migration (Project 1) | `migration_{no_drug,with_drug}_{low,med,high}_density.pkl` | 41 | 5 min | 5 / 48 / 480 | 0 empty, 1 cell |
| infection (Project 2) | `infection_{no_drug,with_drug}_experiment_{0,1,2}.pkl` | 401 | 0.5 min | ~161 | 0 empty, 1 S, 2 I, 3 R |

Field of view is a 40x40 lattice at 1 micron per site. Migration movies are binary
occupancy grids **with no cell identities** — trajectories are reconstructed by
linking (see `src/mism/tracking.py`).

If `data/raw/` is empty after cloning, copy the six `migration_*` (and, for Project 2,
`infection_*`) `.pkl` files into it.
