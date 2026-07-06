# =============================================================================
# tracking.R  --  turn occupancy grids into trajectories, one frame pair at a time.
# =============================================================================
#
# The pipeline is deliberately tiny and readable:
#
#   detect_cells()  a frame  -> the (x, y) coordinates of the occupied sites
#   cost_matrix()   two point sets -> pairwise minimum-image distances
#   build_tracks()  a whole movie -> trajectories, by chaining solve_assignment()
#   diagnostics()   trajectories -> the confidence ratio rho and its ingredients
#   plot_tracks()   trajectories -> a picture
#
# The only interesting idea is in build_tracks(): a "trajectory" is not in the
# data -- it is manufactured by solving one assignment per frame pair and then
# threading each cell's matched index through time. Its reliability is exactly
# the confidence ratio rho = mean step / nearest-neighbour spacing.
#
# Dependencies: base R only (the solver is in hungarian.R). No trajectory magic,
# no black boxes.
# -----------------------------------------------------------------------------

# --- frame access ------------------------------------------------------------
# reticulate may hand back A_list either as a 3-D array (t, x, y) OR as a list of
# 2-D frames, depending on how the pickle was written. These helpers make the
# rest of the code agnostic to that (mirrors extract_frame() in data_share.Rmd).
is_frame_list <- function(A) is.list(A) && !is.array(A)
n_frames <- function(A) if (is_frame_list(A)) length(A) else dim(A)[1]
get_frame <- function(A, t) if (is_frame_list(A)) as.matrix(A[[t]]) else as.matrix(A[t, , ])

# --- detection ---------------------------------------------------------------
# A frame is an L x L integer matrix (0 = empty). Return an M x 2 matrix of the
# occupied-site coordinates. `state = NULL` means "any cell" (value > 0); pass a
# specific code (e.g. 1) to select one cell state.
detect_cells <- function(frame, state = NULL) {
  mask <- if (is.null(state)) frame > 0 else frame == state
  idx  <- which(mask, arr.ind = TRUE)      # columns: row index (x), col index (y)
  colnames(idx) <- NULL
  matrix(as.double(idx), ncol = 2)         # (x, y); absolute origin is irrelevant
}

# --- minimum-image pairwise distance ----------------------------------------
# Distance from every point in `a` (na x 2) to every point in `b` (nb x 2).
# If periodic, wrap each coordinate difference into (-L/2, L/2] so a cell that
# steps across an edge is not mistaken for a giant jump.
cost_matrix <- function(a, b, L, periodic = TRUE) {
  D2 <- matrix(0, nrow(a), nrow(b))
  for (k in 1:2) {
    d <- outer(a[, k], b[, k], "-")        # na x nb difference in dimension k
    if (periodic) d <- ((d + L / 2) %% L) - L / 2
    D2 <- D2 + d^2
  }
  sqrt(D2)
}

# --- link a whole movie into trajectories ------------------------------------
# `A` is the (t, x, y) occupancy array (as returned by reticulate). Returns a
# list with:
#   pos  : (T x N x 2) continuous, unwrapped positions in track order
#   disp : (T-1 x N x 2) one-step minimum-image displacements
# `solver(cost)` must return the column assigned to each row (defaults to our
# from-scratch Hungarian; swap in clue::solve_LSAP for large N).
build_tracks <- function(A, L, periodic = TRUE, solver = solve_assignment, state = NULL) {
  Tn <- n_frames(A)
  p0 <- detect_cells(get_frame(A, 1), state)
  N  <- nrow(p0)

  pos  <- array(NA_real_, dim = c(Tn, N, 2))
  disp <- array(NA_real_, dim = c(Tn - 1, N, 2))
  pos[1, , ] <- p0

  for (t in seq_len(Tn - 1)) {
    cur <- matrix(pos[t, , ], N, 2)              # current positions (may be unwrapped)
    nxt <- detect_cells(get_frame(A, t + 1), state)  # detections next frame

    # cost = distance from each current cell to each detection (wrap cur into
    # the box first; %% keeps the minimum-image cost correct even if unwrapped).
    cost  <- cost_matrix(cur %% L, nxt, L, periodic)
    match <- solver(cost)                        # column (detection) per row (cell)

    for (i in seq_len(N)) {
      d <- nxt[match[i], ] - cur[i, ]            # raw displacement to matched detection
      if (periodic) d <- ((d + L / 2) %% L) - L / 2
      disp[t, i, ]     <- d
      pos[t + 1, i, ]  <- cur[i, ] + d           # keep the trajectory continuous
    }
  }
  list(pos = pos, disp = disp)
}

# --- confidence diagnostics --------------------------------------------------
# rho < ~0.5 : the nearest detection is almost always the true cell (reliable).
# rho ~ 1    : the nearest detection is often the WRONG cell (links are guesses).
diagnostics <- function(tracks, L, periodic = TRUE) {
  disp <- tracks$disp
  step <- sqrt(disp[, , 1]^2 + disp[, , 2]^2)    # (T-1) x N per-step distances
  mean_step <- mean(step, na.rm = TRUE)

  N  <- dim(tracks$pos)[2]
  p0 <- matrix(tracks$pos[1, , ], N, 2)
  Dself <- cost_matrix(p0, p0, L, periodic)
  diag(Dself) <- Inf
  nn <- stats::median(apply(Dself, 1, min))      # median nearest-neighbour spacing

  list(n_cells = N, mean_step = mean_step, nn_spacing = nn, rho = mean_step / nn)
}

# --- draw the reconstructed trajectories -------------------------------------
plot_tracks <- function(tracks, L, main = "", file = NULL) {
  pos <- tracks$pos
  Tn  <- dim(pos)[1]; N <- dim(pos)[2]
  xr  <- range(pos[, , 1]); yr <- range(pos[, , 2])
  if (!is.null(file)) png(file, width = 900, height = 900, res = 130)
  plot(NA, xlim = xr, ylim = yr, xlab = "x (microns)", ylab = "y (microns)",
       main = main, asp = 1)
  cols <- grDevices::rainbow(N)
  for (i in seq_len(N)) {
    lines(pos[, i, 1], pos[, i, 2], col = cols[i])
    points(pos[1, i, 1],  pos[1, i, 2],  col = cols[i], pch = 1, cex = 1.1)   # start
    points(pos[Tn, i, 1], pos[Tn, i, 2], col = cols[i], pch = 16, cex = 1.0)  # end
  }
  if (!is.null(file)) { dev.off(); message("wrote ", file) }
}

# --- tidy exports ------------------------------------------------------------
# The trajectories as a long ("tidy") data frame: one row per (cell, frame).
# This is `tr$pos` unrolled -- cell i is the same cell at every frame.
#
# Pass the raw movie `A` to also attach each cell's STATE. Tracking only assigns
# identity (the cell index); the S/I/R value is read verbatim from the frames at
# the cell's tracked lattice site -- it is ground truth from the data, not
# inferred. (For the 0/1 migration movies the state is just 1 = cell.)
tracks_to_df <- function(tracks, dt = 1, A = NULL, L = NULL) {
  pos <- tracks$pos
  Tn <- dim(pos)[1]; N <- dim(pos)[2]
  df <- data.frame(
    cell  = rep(seq_len(N), each = Tn),
    frame = rep(seq_len(Tn), times = N),
    t_min = rep(seq_len(Tn) - 1, times = N) * dt,
    x     = as.vector(pos[, , 1]),   # column-major: frame varies fastest within a cell
    y     = as.vector(pos[, , 2])
  )
  if (!is.null(A)) {
    if (is.null(L)) L <- nrow(get_frame(A, 1))
    # tracked positions are integer lattice sites (possibly unwrapped); wrap them
    # back into 1..L and read the state at fr[row, col] = fr[x, y].
    xi <- ((round(df$x) - 1) %% L) + 1
    yi <- ((round(df$y) - 1) %% L) + 1
    st <- integer(nrow(df))
    for (t in seq_len(Tn)) {
      fr  <- get_frame(A, t)
      sel <- which(df$frame == t)
      st[sel] <- fr[cbind(xi[sel], yi[sel])]
    }
    df$state       <- st
    df$state_label <- factor(st, levels = 0:3, labels = c("empty", "S", "I", "R"))
  }
  df
}

# The LINKS as a long data frame: one row per (cell, frame t -> t+1). Each row is
# a single frame-to-frame match, so filtering to one `frame_from` gives you the
# links for exactly one step.
links_to_df <- function(tracks, dt = 1) {
  pos <- tracks$pos; disp <- tracks$disp
  K <- dim(disp)[1]; N <- dim(disp)[2]          # K = number of steps = T - 1
  x0 <- as.vector(pos[seq_len(K), , 1]); y0 <- as.vector(pos[seq_len(K), , 2])
  dx <- as.vector(disp[, , 1]);          dy <- as.vector(disp[, , 2])
  data.frame(
    cell       = rep(seq_len(N), each = K),
    frame_from = rep(seq_len(K), times = N),
    frame_to   = rep(seq_len(K), times = N) + 1L,
    t_min      = rep(seq_len(K) - 1, times = N) * dt,
    x0 = x0, y0 = y0, x1 = x0 + dx, y1 = y0 + dy,
    dx = dx, dy = dy, dist = sqrt(dx^2 + dy^2)
  )
}

# --- inspect a SINGLE link (one frame pair) ----------------------------------
# Draw only the matching for frame t -> t+1: open circles at t, filled dots at
# t+1, and an arrow along each cell that moved. Contrast with plot_tracks(),
# which overlays the whole multi-frame trajectory.
plot_link <- function(tracks, t, cells = NULL, main = NULL, file = NULL, pad = 1) {
  pos <- tracks$pos
  a <- matrix(pos[t, , ],     ncol = 2)
  b <- matrix(pos[t + 1, , ], ncol = 2)
  if (is.null(cells)) cells <- seq_len(nrow(a))
  a <- a[cells, , drop = FALSE]; b <- b[cells, , drop = FALSE]
  moved <- rowSums((b - a)^2) > 0
  if (is.null(main))
    main <- sprintf("Link  frame %d -> %d   (%d cells, %d moved)",
                    t, t + 1, nrow(a), sum(moved))
  xr <- range(c(a[, 1], b[, 1])) + c(-pad, pad)
  yr <- range(c(a[, 2], b[, 2])) + c(-pad, pad)
  if (!is.null(file)) png(file, width = 900, height = 900, res = 130)
  plot(NA, xlim = xr, ylim = yr, xlab = "x (microns)", ylab = "y (microns)",
       main = main, asp = 1)
  points(a[, 1], a[, 2], pch = 1,  col = "grey45",   cex = 1.2)   # frame t
  points(b[, 1], b[, 2], pch = 16, col = "steelblue", cex = 0.9)  # frame t+1
  if (any(moved))
    suppressWarnings(arrows(a[moved, 1], a[moved, 2], b[moved, 1], b[moved, 2],
                            length = 0.07, col = "firebrick", lwd = 1.5))
  legend("topright", c(sprintf("frame %d", t), sprintf("frame %d", t + 1), "link"),
         pch = c(1, 16, NA), lty = c(NA, NA, 1),
         col = c("grey45", "steelblue", "firebrick"), bty = "n")
  if (!is.null(file)) { dev.off(); message("wrote ", file) }
}
