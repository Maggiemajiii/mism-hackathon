# MISM Summer 2026 Hackathon

Organized workspace for the MISM Summer 2026 hackathon: quantifying **cell movement**
(Project 1) and **antiviral effect on infection spread** (Project 2) from time-lapse
microscopy stored as occupancy grids.

## Layout

```
data/
  raw/                 pickled movies (gitignored; kept locally)
  README.md            data dictionary (format + cell-state encoding)
docs/                  assignment brief (MISM_hackathon_summer2026.pdf)
notebooks/             organizers' starter notebooks (data_share.ipynb / .Rmd / .pdf)
src/                   analysis code (to be added)
project/               per-project workspaces (to be added)
MISM_workshop.Rproj    RStudio project file
```

## Setup

```bash
# Python
pip install -r requirements.txt        # or: pip install -e .
# R
# open MISM_workshop.Rproj in RStudio
```

## Data

Twelve movies live in `data/raw/` (see `data/README.md` for the full dictionary):

- `migration_{no_drug,with_drug}_{low,med,high}_density.pkl` — 41 frames, 5 min interval
- `infection_{no_drug,with_drug}_experiment_{0,1,2}.pkl` — 401 frames, 0.5 min interval

Each is a 40x40 lattice at 1 micron/site. `.pkl` files are gitignored (kept locally,
not committed).
