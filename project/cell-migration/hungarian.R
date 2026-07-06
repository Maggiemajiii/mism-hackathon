# =============================================================================
# hungarian.R  --  the assignment problem, spelled out from scratch.
# =============================================================================
#
# WHY THIS EXISTS
# ---------------
# Our movies are occupancy grids: at frame t we know *where* the cells are, but
# not *which* cell is which. To turn that into trajectories we have to link each
# cell at frame t to its future self at frame t+1. Phrase it as a matching:
#
#     given N cells (rows) and their N detections next frame (columns),
#     pick a one-to-one matching row -> column that minimizes total movement.
#
# Formally, choose a permutation `match` to
#
#     minimize   sum_i  cost[i, match(i)]
#     s.t.       each row and each column is used exactly once.
#
# That is the classic *linear assignment problem*. The Hungarian algorithm
# (Kuhn, 1955) solves it EXACTLY in O(n^3). It is a primal-dual method: it keeps
# dual "potentials" u_i (rows) and v_j (columns), and on each outer step it grows
# a shortest alternating path and pulls one more row into the matching. Greedy
# nearest-neighbour can double-book a detection; this cannot -- it is optimal.
#
# WHAT YOU GET
# ------------
#   solve_assignment(cost)  ->  integer vector `match`, length nrow(cost),
#   where match[i] is the column assigned to row i. For a square cost matrix
#   (our case: #cells == #detections) it is a permutation of 1..n.
#
# This 1-indexed "potentials + shortest augmenting path" version was checked
# against scipy.optimize.linear_sum_assignment on 3000 random matrices: the
# total assignment cost matched every time.
#
# It is O(n^3) in plain R loops -- perfect for the low/medium-density movies
# (5 and 48 cells). For the 480-cell movie prefer clue::solve_LSAP (same
# problem, C speed); run_baseline.R falls back to it automatically.
# -----------------------------------------------------------------------------

solve_assignment <- function(cost) {
  n <- nrow(cost)          # rows    = cells at frame t
  m <- ncol(cost)          # columns = detections at frame t+1  (n <= m)
  stopifnot(n <= m)

  # --- 1-indexed bookkeeping ------------------------------------------------
  # There is a virtual "sentinel" column 0 that temporarily holds the row being
  # added. We store logical column c (0..m) at vector position c+1, so position
  # 1 is the sentinel. Likewise u[r+1] is the potential of logical row r.
  u   <- numeric(n + 1)    # row potentials
  v   <- numeric(m + 1)    # column potentials
  p   <- integer(m + 1)    # p[c+1] = row currently matched to column c (0 = free)
  way <- integer(m + 1)    # back-pointers describing the alternating path

  for (i in seq_len(n)) {          # add row i to the matching
    p[1] <- i                      # park it on the sentinel column (logical 0)
    j0 <- 0L                       # current column of the search
    minv <- rep(Inf, m + 1)        # cheapest reduced cost reaching each column
    used <- rep(FALSE, m + 1)

    repeat {                       # Dijkstra-like growth of the alternating tree
      used[j0 + 1] <- TRUE
      i0    <- p[j0 + 1]           # the row sitting on column j0
      delta <- Inf
      j1    <- -1L
      for (j in seq_len(m)) {      # relax every not-yet-used real column
        if (!used[j + 1]) {
          cur <- cost[i0, j] - u[i0 + 1] - v[j + 1]   # reduced cost
          if (cur < minv[j + 1]) { minv[j + 1] <- cur; way[j + 1] <- j0 }
          if (minv[j + 1] < delta) { delta <- minv[j + 1]; j1 <- j }
        }
      }
      for (j in 0:m) {             # shift potentials by delta (keeps duals valid)
        if (used[j + 1]) {
          u[p[j + 1] + 1] <- u[p[j + 1] + 1] + delta
          v[j + 1]        <- v[j + 1] - delta
        } else {
          minv[j + 1] <- minv[j + 1] - delta
        }
      }
      j0 <- j1
      if (p[j0 + 1] == 0) break    # reached a free column -> we can augment
    }

    repeat {                       # walk the back-pointers, flipping the matches
      j1 <- way[j0 + 1]
      p[j0 + 1] <- p[j1 + 1]
      j0 <- j1
      if (j0 == 0) break
    }
  }

  # p[c+1] holds the row matched to column c; invert it to "column per row".
  match <- integer(n)
  for (j in seq_len(m)) if (p[j + 1] != 0) match[p[j + 1]] <- j
  match
}
