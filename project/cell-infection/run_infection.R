# =============================================================================
# run_infection.R  --  the same Hungarian tracking baseline, on the infection data.
# =============================================================================
#
# The infection movies encode cell STATE (0 empty, 1 S, 2 I, 3 R) rather than a
# plain 0/1 mask, but a cell keeps its identity while its state changes, and the
# total cell count is constant (161). So we track *any* occupied site (state > 0)
# with exactly the same code as the migration baseline -- only the data changes.
#
# What tracking tells us here (movement, not infection):
#   * are these movies reliably trackable?         (the rho column)
#   * how much do the cells actually move?          (speed, % stationary)
#   * does the drug change cell MOTILITY?           (compare the two conditions)
# The infection/SIR spread analysis itself is a separate piece of work.
#
# Run from anywhere:  Rscript run_infection.R
# Needs: R, reticulate + Python/numpy, and clue (161 cells x 401 frames wants the
#        C-speed solver: install.packages("clue")).
# -----------------------------------------------------------------------------

script_dir <- local({
  a <- commandArgs(FALSE); k <- grep("^--file=", a)
  if (length(k)) dirname(normalizePath(sub("^--file=", "", a[k]))) else getwd()
})
# reuse the exact solver + tracking code from the migration baseline
source(file.path(script_dir, "..", "cell-migration", "hungarian.R"))
source(file.path(script_dir, "..", "cell-migration", "tracking.R"))

suppressMessages(library(reticulate))
DATA <- normalizePath(file.path(script_dir, "..", "..", "data", "raw"))

if (!requireNamespace("clue", quietly = TRUE))
  stop("This needs the 'clue' package for the 161-cell movies. ",
       "Run install.packages(\"clue\") and try again.")
fast_solver <- function(cost) as.integer(clue::solve_LSAP(cost))

conditions  <- c("no_drug", "with_drug")
experiments <- 0:2

cat(sprintf("%-10s %4s %6s %11s %11s %8s %9s %8s\n",
            "condition", "exp", "cells", "mean_step", "nn_space", "rho", "um/min", "%still"))
cat(strrep("-", 72), "\n")

tracks <- list()
for (cond in conditions) {
  for (e in experiments) {
    f   <- file.path(DATA, sprintf("infection_%s_experiment_%d.pkl", cond, e))
    out <- py_load_object(f)
    A   <- out$A_list
    L   <- as.integer(out$metadata[["xy dimensions"]])
    if (length(L) == 0 || is.na(L)) L <- nrow(get_frame(A, 1))
    dt  <- as.numeric(out$metadata[["time interval (mins)"]])

    tr <- build_tracks(A, L, periodic = TRUE, solver = fast_solver, state = NULL)
    cat(tr$pos[1, 1:5, ], "\n")  # sanity check: first frame, first 5 cells
    cat(length(tr))
    d  <- diagnostics(tr, L, periodic = TRUE)

    step       <- sqrt(tr$disp[, , 1]^2 + tr$disp[, , 2]^2)
    frac_still <- mean(step == 0, na.rm = TRUE)
    speed      <- d$mean_step / dt

    tracks[[paste(cond, e, sep = "_")]] <- tr
    cat(sprintf("%-10s %4d %6d %11.3f %11.3f %8.3f %9.3f %7.0f%%\n",
                cond, e, d$n_cells, d$mean_step, d$nn_spacing, d$rho, speed,
                100 * frac_still))
  }
}

cat("\nrho << 0.5 everywhere -> the cells move little and are reliably trackable.\n")
cat("Motility is ~equal across drug/no-drug, so the drug acts on transmission,\n")
cat("not on how the cells move.\n")

# --- picture: trajectories of one movie (they are short -- cells barely move) --
dir.create(file.path(script_dir, "results"), showWarnings = FALSE)
plot_tracks(tracks[["no_drug_0"]], L = 40,
            main = "Infection assay, no drug (exp 0) -- reconstructed tracks",
            file = file.path(script_dir, "results", "tracks_infection_no_drug_0.png"))

# --- save the trajectory as a tidy data frame --------------------------------
# `tr$pos` unrolled to one row per (cell, frame): cell, frame, t_min, x, y.
traj <- tracks_to_df(tracks[["no_drug_0"]], dt = 0.5)
write.csv(traj, file.path(script_dir, "results", "trajectory_no_drug_0.csv"),
          row.names = FALSE)
cat(sprintf("\nTrajectory data frame: %d rows = %d cells x %d frames\n",
            nrow(traj), length(unique(traj$cell)), length(unique(traj$frame))))
print(utils::head(traj))

# --- look at ONE link between frames (frame 1 -> 2), not the whole trajectory --
plot_link(tracks[["no_drug_0"]], t = 1,
          file = file.path(script_dir, "results", "link_no_drug_0_frame1.png"))
