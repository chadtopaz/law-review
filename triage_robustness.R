################################################################################
# Triage Robustness Check
#
# Requested by Reviewer 1 (Major Weakness 1) on the PLOS ONE revision of
# PONE-D-26-11974. The reviewer notes that real student editors likely
# desk-reject many submissions on a cursory glance, reserving deep review
# capacity for promising manuscripts. The current model captures this only
# partially through quality-dependent review effort (Eq. 3 / compute_review_count);
# it lacks an explicit desk-rejection step.
#
# This script implements an explicit triage step and reports a focused
# Monte Carlo comparison: Current System vs. Current System + Triage vs.
# Deferred Acceptance. All other mechanism details are unchanged.
#
# Triage design:
#   - Triggered once per journal at the start of the simulation (since all
#     articles arrive on Day 1 under batch submission).
#   - Each journal screens its entire initial queue at negligible cost that
#     does not consume the journal's regular review capacity c_j. The
#     screening uses the same signal q_hat_ji that drives deep review.
#   - Articles whose signal falls below the journal's TRIAGE_THRESHOLD_PCT
#     percentile of its queue are removed from that journal's queue
#     (desk-rejected at THIS journal). They remain in other journals' queues.
#   - This follows the referee's framing: triage is a cheap screening step
#     and c_j is reserved for the surviving pool. We model the cheap step
#     as zero-binding-cost rather than as an explicit cost that consumes a
#     separate budget; under realistic capacity multipliers, an explicit
#     cost would not bind in any case, so the simplification is innocuous.
#   - The main simulation loop then proceeds identically to simulate_current.
#
# Output: output/triage_robustness.rds
################################################################################

# Load core simulation (do not auto-run the full pipeline).
# Look in both the full-project layout (code/simulation.R) and the
# public-replication-repo layout (simulation.R at root) so this script
# works whether run from a clone of the public repo or from the author's
# full working project.
#
# We save and restore the caller's RUN_PIPELINE value to avoid accidentally
# triggering the full baseline pipeline if the user previously set
# RUN_PIPELINE <- TRUE in the same R session.
.prev_run_pipeline <- if (exists("RUN_PIPELINE")) RUN_PIPELINE else NULL
RUN_PIPELINE <- FALSE
if (file.exists("simulation.R")) {
  source("simulation.R")
} else if (file.exists("code/simulation.R")) {
  source("code/simulation.R")
} else {
  stop("Cannot find simulation.R. Run from the project root or from a ",
       "clone of github.com/chadtopaz/law-review.")
}
if (is.null(.prev_run_pipeline)) {
  rm(RUN_PIPELINE)
} else {
  RUN_PIPELINE <- .prev_run_pipeline
}
rm(.prev_run_pipeline)


# ==============================================================================
# simulate_current_triage(): variant of simulate_current with explicit triage
# ==============================================================================

simulate_current_triage <- function(journals, articles, p,
                                     signal_mat = NULL,
                                     crn = NULL,
                                     triage_threshold_pct = 0.40) {
  J <- p$J
  N <- p$N
  T_steps <- p$T_steps
  deadline_accept <- identical(p$deadline_behavior, "accept")

  if (is.null(signal_mat)) signal_mat <- make_signal_matrix(articles, p, journals)
  if (is.null(crn)) crn <- make_crn_draws(articles, p)
  j_queue <- lapply(crn$j_queue, function(x) x)  # deep copy queues

  # ============================================================================
  # TRIAGE STEP (pre-simulation)
  # Each journal screens its entire queue. Articles below the journal's
  # threshold percentile of queue signals are desk-rejected (removed from
  # this journal's queue only).
  # ============================================================================
  for (jj in 1:J) {
    q <- j_queue[[jj]]
    if (length(q) == 0) next
    sigs <- signal_mat[jj, q]
    # Journal-specific threshold: percentile of queue signal
    thr <- quantile(sigs, probs = triage_threshold_pct, names = FALSE)
    survivors <- q[sigs >= thr]
    j_queue[[jj]] <- survivors
  }

  # ============================================================================
  # Main loop: identical to simulate_current
  # ============================================================================
  j_remaining     <- journals$slots
  a_best_offer    <- rep(Inf, N)
  a_best_j        <- rep(NA_integer_, N)
  a_matched_j     <- rep(NA_integer_, N)
  a_status        <- rep(STATUS_SUBMITTED, N)
  a_offer_expires <- rep(Inf, N)

  for (t in 1:T_steps) {

    # Handle expired offers
    expired <- which(a_offer_expires <= t & a_status == STATUS_HOLDING)
    if (length(expired) > 0) {
      if (deadline_accept) {
        a_matched_j[expired] <- a_best_j[expired]
        a_status[expired] <- STATUS_MATCHED
      } else {
        for (ei in expired) {
          j_remaining[a_best_j[ei]] <- j_remaining[a_best_j[ei]] + 1L
          a_best_offer[ei] <- Inf
          a_best_j[ei] <- NA_integer_
          a_status[ei] <- STATUS_SUBMITTED
          a_offer_expires[ei] <- Inf
        }
      }
    }

    j_order <- crn$j_orders[t, ]

    for (jj in j_order) {
      j_committed <- journals$slots[jj] - j_remaining[jj]
      overbook_cap <- ceiling(journals$slots[jj] * journals$overbook[jj])
      if (j_committed >= overbook_cap) next
      if (length(j_queue[[jj]]) == 0) next

      queue <- j_queue[[jj]]
      queue <- queue[a_status[queue] != STATUS_MATCHED]
      j_queue[[jj]] <- queue
      if (length(queue) == 0) next

      # Expedite prioritization
      is_expedite <- a_status[queue] == STATUS_HOLDING
      if (any(is_expedite) && crn$expedite_coins[t, jj] <= journals$expedite_prob[jj]) {
        queue_sorted <- c(queue[is_expedite], queue[!is_expedite])
      } else {
        queue_sorted <- queue
      }

      # Quality-dependent review with c_j (unchanged from simulate_current)
      n_review <- compute_review_count(queue_sorted, journals$capacity[jj],
                                        signal_mat, jj)
      n_review <- min(n_review, length(queue_sorted))
      reviewed <- queue_sorted[1:n_review]
      j_queue[[jj]] <- queue_sorted[-(1:n_review)]

      sigs <- signal_mat[jj, reviewed]

      n_offers <- min(overbook_cap - j_committed, length(reviewed))
      n_offers <- max(n_offers, 0L)
      if (n_offers == 0) next

      top_idx <- order(sigs, decreasing = TRUE)[1:n_offers]
      offer_articles <- reviewed[top_idx]

      for (ai in offer_articles) {
        if (a_status[ai] == STATUS_MATCHED) next

        if (jj < a_best_offer[ai]) {
          old_j <- a_best_j[ai]

          a_best_offer[ai] <- jj
          a_best_j[ai] <- jj
          a_status[ai] <- STATUS_HOLDING
          a_offer_expires[ai] <- t + journals$deadline[jj]

          if (!is.na(old_j)) {
            j_remaining[old_j] <- j_remaining[old_j] + 1L
          }
          j_remaining[jj] <- j_remaining[jj] - 1L

          time_frac <- t / T_steps
          accept_pct <- 0.10 + 0.50 * time_frac
          top_threshold <- articles$aspiration[ai] +
            accept_pct * (articles$sub_hi[ai] - articles$aspiration[ai])
          if (jj <= floor(top_threshold)) {
            a_matched_j[ai] <- jj
            a_status[ai] <- STATUS_MATCHED
          }
        }
      }
    }
  }

  still_holding <- which(a_status == STATUS_HOLDING)
  a_matched_j[still_holding] <- a_best_j[still_holding]

  a_matched_j <- resolve_overbooking(a_matched_j, journals, articles, signal_mat, J)

  data.table(
    a_id      = 1:N,
    quality   = articles$quality,
    prestige  = articles$prestige,
    matched_j = a_matched_j
  )
}


# ==============================================================================
# run_triage_one_rep(): per-replication work (top-level for cluster serialization)
# Following the same global-helper + inline-wrapper pattern as simulation.R's
# run_monte_carlo() to avoid closure-serialization issues across workers.
# ==============================================================================

run_triage_one_rep <- function(i, p, triage_threshold_pct) {
  journals <- init_journals(p)
  articles <- init_articles(p)
  signal_mat <- make_signal_matrix(articles, p, journals)
  crn <- make_crn_draws(articles, p)

  res_current <- simulate_current(journals, articles, p,
                                   signal_mat = signal_mat, crn = crn)
  res_triage  <- simulate_current_triage(journals, articles, p,
                                          signal_mat = signal_mat, crn = crn,
                                          triage_threshold_pct = triage_threshold_pct)
  res_da      <- simulate_da(journals, articles, p, proposer = "journal",
                              signal_mat = signal_mat)

  list(
    current = compute_welfare(res_current, journals, articles, p),
    triage  = compute_welfare(res_triage,  journals, articles, p),
    da      = compute_welfare(res_da,      journals, articles, p)
  )
}


# ==============================================================================
# run_triage_robustness(): focused Monte Carlo
# ==============================================================================

run_triage_robustness <- function(p = params, n_sims = 500, ncores = NCORES,
                                   triage_threshold_pct = 0.40) {
  # Force evaluation of all arguments BEFORE any parallel work.
  # Without this, the default value `p = params` is a lazy promise that
  # would be evaluated only when first read -- and if that read happens on
  # a worker node (which lacks `params`), it throws "object 'params' not found".
  force(p)
  force(n_sims)
  force(ncores)
  force(triage_threshold_pct)

  message("=========================================================")
  message("Triage Robustness Check (Reviewer 1 Major Weakness 1)")
  message("=========================================================")
  message("  n_sims:                ", n_sims)
  message("  triage_threshold_pct:  ", triage_threshold_pct,
          "  (reject signals below the ", round(100*triage_threshold_pct),
          "th percentile of each journal's queue;",
          " screening is negligible-cost and does not consume c_j)")
  message("")

  cl <- makeCluster(ncores)
  on.exit(stopCluster(cl))

  clusterExport(cl, c(
    "init_journals", "init_articles",
    "make_signal_matrix", "build_queues", "make_crn_draws",
    "resolve_overbooking", "compute_review_count",
    "simulate_current", "simulate_da", "simulate_current_triage",
    "compute_welfare", "run_triage_one_rep",
    "STATUS_SUBMITTED", "STATUS_HOLDING", "STATUS_MATCHED"
  ))
  clusterEvalQ(cl, library(data.table))
  clusterSetRNGStream(cl, iseed = 2026)

  t0 <- Sys.time()
  # Inline anonymous wrapper captures p and triage_threshold_pct from
  # this scope and passes them explicitly to the global helper.
  # This matches the pattern used in simulation.R's run_monte_carlo().
  results <- pblapply(1:n_sims, function(i) {
    run_triage_one_rep(i, p, triage_threshold_pct)
  }, cl = cl)
  t1 <- Sys.time()
  message("\nMonte Carlo completed in ",
          round(difftime(t1, t0, units = "mins"), 1), " minutes.")

  # ----- Aggregate -----
  mechanisms <- c("current", "triage", "da")
  mech_labels <- c("Current System", "Current + Triage", "Deferred Acceptance")

  agg <- data.table(
    mechanism      = mech_labels,
    match_quality  = sapply(mechanisms, function(m)
      mean(sapply(results, function(r) r[[m]]$match_quality))),
    rank_cor       = sapply(mechanisms, function(m)
      mean(sapply(results, function(r) r[[m]]$rank_cor))),
    n_matched      = sapply(mechanisms, function(m)
      mean(sapply(results, function(r) r[[m]]$n_matched))),
    all_N_welfare  = sapply(mechanisms, function(m)
      mean(sapply(results, function(r) r[[m]]$match_quality_allN))),
    tier1_quality  = sapply(mechanisms, function(m)
      mean(sapply(results, function(r) r[[m]]$tier_avg_article_q[1]))),
    tier2_quality  = sapply(mechanisms, function(m)
      mean(sapply(results, function(r) r[[m]]$tier_avg_article_q[2]))),
    tier3_quality  = sapply(mechanisms, function(m)
      mean(sapply(results, function(r) r[[m]]$tier_avg_article_q[3]))),
    tier4_quality  = sapply(mechanisms, function(m)
      mean(sapply(results, function(r) r[[m]]$tier_avg_article_q[4])))
  )

  # Gap vs. DA on all-N welfare (positive = DA better)
  da_q <- agg[mechanism == "Deferred Acceptance"]$all_N_welfare
  agg[, gap_vs_da_pct := round(100 * (all_N_welfare - da_q) / da_q, 1)]

  message("\n--- Triage Robustness Results ---")
  print(agg)

  # ----- Save -----
  out <- list(
    aggregate            = agg,
    raw                  = results,
    triage_threshold_pct = triage_threshold_pct,
    n_sims               = n_sims,
    elapsed_minutes      = as.numeric(difftime(t1, t0, units = "mins"))
  )

  if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
  saveRDS(out, file.path(OUTPUT_DIR, "triage_robustness.rds"))
  message("Saved to ", file.path(OUTPUT_DIR, "triage_robustness.rds"))

  invisible(out)
}


# ==============================================================================
# Execute (set SKIP_RUN <- TRUE before sourcing to suppress)
# ==============================================================================
if (!exists("SKIP_RUN") || !isTRUE(SKIP_RUN)) {
  triage_results <- run_triage_robustness()
}
