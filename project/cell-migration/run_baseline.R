# =============================================================================
# run_baseline.R  --  baseline cell-trajectory reconstruction (Hungarian linking)
# =============================================================================
#
# Loads each migration movie, reconstructs trajectories by frame-to-frame
# Hungarian assignment, prints the confidence table, and draws the trajectories
# for the one movie we can actually trust (low density, no drug).
#
# Run from anywhere:   Rscript run_baseline.R
# Needs: R, reticulate + a Python with numpy (same setup as data_share.Rmd).
#        clue is optional -- only used to speed up the 480-cell movie.
# -----------------------------------------------------------------------------

# --- locate ourselves so sourcing + data paths work from any working dir -----
script_dir <- local({
  a <- commandArgs(FALSE); k <- grep("^--file=", a)
  if (length(k)) dirname(normalizePath(sub("^--file=", "", a[k]))) else getwd()
})
source(file.path(script_dir, "hungarian.R"))
source(file.path(script_dir, "tracking.R"))

suppressMessages(library(reticulate))
DATA <- normalizePath(file.path(script_dir, "..", "..", "data", "raw"))

# A fast solver for big movies: clue::solve_LSAP solves the same assignment in C.
fast_solver <- if (requireNamespace("clue", quietly = TRUE)) {
  function(cost) as.integer(clue::solve_LSAP(cost))
} else NULL

conditions <- c("no_drug", "with_drug")
densities  <- c("low", "med", "high")

cat(sprintf("%-10s %-5s %6s %11s %11s %8s\n",
            "condition", "dens", "cells", "mean_step", "nn_space", "rho"))
cat(strrep("-", 54), "\n")

tracks <- list()
for (cond in conditions) {
  for (dens in densities) {
    f <- file.path(DATA, sprintf("migration_%s_%s_density.pkl", cond, dens))
    out <- py_load_object(f)                       # dict -> named list
    A   <- out$A_list                              # 3-D array OR list of frames
    L   <- as.integer(out$metadata[["xy dimensions"]])
    if (length(L) == 0 || is.na(L)) L <- nrow(get_frame(A, 1))

    n_cells <- nrow(detect_cells(get_frame(A, 1)))

    # from-scratch Hungarian is O(n^3): fine for small N. For the 480-cell movie
    # use clue::solve_LSAP (C speed); if clue is absent, skip it rather than hang
    # (run install.packages("clue") to include the high-density rows).
    if (n_cells > 150 && is.null(fast_solver)) {
      cat(sprintf("%-10s %-5s %6d %11s %11s %8s\n",
                  cond, dens, n_cells, "-", "-", "skip(need clue)"))
      next
    }
    solver <- if (n_cells > 150) fast_solver else solve_assignment

    tr <- build_tracks(A, L, periodic = TRUE, solver = solver)
    d  <- diagnostics(tr, L, periodic = TRUE)
    tracks[[paste(cond, dens, sep = "_")]] <- tr

    tag <- if (n_cells > 150) "  [clue]" else ""
    cat(sprintf("%-10s %-5s %6d %11.3f %11.3f %8.3f%s\n",
                cond, dens, d$n_cells, d$mean_step, d$nn_spacing, d$rho, tag))
  }
}

cat("\nrho = mean per-frame step / median nearest-neighbour spacing.\n")
cat("rho < ~0.5 : trajectories are reliable;  rho ~ 1 : links are guesses.\n")

# --- picture: the reliably tracked movie -------------------------------------
dir.create(file.path(script_dir, "results"), showWarnings = FALSE)
plot_tracks(tracks[["no_drug_low"]], L = 40,
            main = "Baseline Hungarian tracks -- no drug, low density",
            file = file.path(script_dir, "results", "tracks_no_drug_low.png"))
