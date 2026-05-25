################################################################################
# Law Review Submission System: Agent-Based Simulation
#
# Research Questions:
#   1. What is the magnitude of match quality loss produced by the current
#      decentralized expedite system vs. plausible alternatives?
#   2. How does this inefficiency distribute across the prestige hierarchy?
#
# Author: Chad M. Topaz
# Date: March 2026
################################################################################

# ==============================================================================
# 0. SETUP AND DEPENDENCIES
# ==============================================================================

required_packages <- c("parallel", "pbapply", "data.table", "ggplot2", "scales", "patchwork")
missing <- required_packages[!required_packages %in% rownames(installed.packages())]
if (length(missing) > 0) {
  message("Installing: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
}

library(parallel)
library(pbapply)
library(data.table)
library(ggplot2)
library(scales)
library(patchwork)

NCORES <- detectCores() - 1
FIGURE_DIR <- file.path("manuscript", "plos", "figures")
OUTPUT_DIR <- "output"
message("Using ", NCORES, " cores for parallel computation.")

# --- Figure styling constants ---
# Colorblind-friendly palette (Okabe-Ito). Used only where color is needed.
CB_PALETTE <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
CB_SHAPES <- c(16, 17, 15, 18, 8, 4, 3)  # distinguishable in B&W

# Two-group palettes for Current System vs DA comparisons
TWO_GROUP_FILL <- scale_fill_manual(values = c(
  "Current System" = "grey35", "Deferred Acceptance" = "grey75"))
TWO_GROUP_COLORS <- c("Current System" = "grey25", "Deferred Acceptance" = "grey60")
TWO_GROUP_SHAPES <- c("Current System" = 16, "Deferred Acceptance" = 17)

# PLOS ONE figure specs: max width 7.5in, max height 8.75in, 300-600 DPI.
# Typeset text column width is 5.2in; we use that for single-column figures.
# Output as TIFF with LZW compression; filenames: Fig1.tif, Fig2.tif, etc.
FIG_W <- 5.2; FIG_H <- 3.6; FIG_DPI <- 300
FIG_FORMAT <- "tiff"
FIG_EXT <- ".tif"

FIG_THEME <- theme_bw(base_size = 11, base_family = "Arial") +
  theme(text = element_text(size = 11),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        plot.title = element_text(size = 11, face = "bold"),
        panel.grid.minor = element_blank())

# Helper: save a figure as PLOS-compliant TIFF
save_plos_fig <- function(filename, plot, width = FIG_W, height = FIG_H, dpi = FIG_DPI) {
  ggsave(filename, plot,
         device = FIG_FORMAT, width = width, height = height,
         units = "in", dpi = dpi, compression = "lzw")
}

set.seed(2026)

# ==============================================================================
# Integer status codes for article state tracking
# ==============================================================================
STATUS_SUBMITTED <- 1L
STATUS_HOLDING   <- 2L
STATUS_MATCHED   <- 3L

# ==============================================================================
# 1. MODEL PARAMETERS
# ==============================================================================

params <- list(
  # --- Population sizes ---
  J           = 200,        # Number of journals
  N           = 800,        # Number of articles per cycle

  # --- Quality distribution ---
  q_alpha     = 2,          # Beta distribution shape1 for article quality
  q_beta      = 5,          # Beta distribution shape2 (right-skewed)

  # --- Epsilon (prestige noise) distribution ---
  eps_alpha   = 2,          # Beta distribution shape1 for epsilon
  eps_beta    = 5,          # Beta distribution shape2 for epsilon

  # --- Prestige-quality relationship ---
  # NOTE: This is a mixing weight, not the Pearson correlation itself.
  # Actual cor depends on Var(quality) and Var(epsilon); see compute_actual_rho().
  # For default Beta(2,5) on both, rho=0.4 => actual correlation ~0.555.
  rho         = 0.4,        # Mixing weight for prestige = rho*q + (1-rho)*eps

  # --- Signal structure ---
  w           = 0.6,        # Fallback signal weight (used only if journals lack w_j;
                            # init_journals() always creates journal-specific w_j)
  sigma       = 0.15,       # Noise sd in editorial signal

  # --- Journal parameters ---
  # Slots by tier (min, max) — per-cycle capacity.
  # Calibrated to match annual publication rates of ~4-20 unsolicited articles
  # (Manley 2017) across two main submission seasons.
  slots_t1    = c(3, 6),    # Tier 1: ranks 1-20  (annual ~6-12)
  slots_t2    = c(4, 8),    # Tier 2: ranks 21-50 (annual ~8-16)
  slots_t3    = c(5, 10),   # Tier 3: ranks 51-100 (annual ~10-20)
  slots_t4    = c(6, 12),   # Tier 4: ranks 101-200 (annual ~12-24)

  # Review capacity per time step
  cap_t1      = c(30, 50),
  cap_t2      = c(20, 40),
  cap_t3      = c(10, 30),
  cap_t4      = c(10, 30),

  # Offer deadlines (in time steps)
  # Deadlines vary sharply by journal and offer context across ALL tiers.
  # Elite journals can use very short deadlines on expedited-review offers
  # (e.g., Penn 1-hour, Columbia 1-hour) or longer windows (HLR 7-day).
  # Lower-ranked journals also use short exploding offers. We use one
  # set for Tier 1 and a different (shorter) set for other tiers;
  # the key variation is within-tier, not simply between tiers.
  deadline_t1 = c(1, 3, 5, 7),
  deadline_other = c(1, 1, 2, 3),

  # Expedite responsiveness: probability that a journal honors an expedite
  # request by moving the article to the front of its queue. Not all
  # journals treat expedites uniformly — some say expedites confer no
  # advantage (e.g., Penn, NYU JILP), others attempt to honor them.
  # Drawn journal-by-journal from Beta(mean * conc, (1-mean) * conc).
  # No systematic tier gradient: evidence does not clearly support
  # "more elite = more responsive" (Penn says no advantage; lower-tier
  # journals may be more eager to attract authors via expedite).
  expedite_prob_mean = 0.65,  # Population mean across all journals
  expedite_prob_conc = 8,     # Concentration (sd ~ 0.16); moderate spread

  # Overbooking factor
  overbook_t1 = 1.5,
  overbook_t2 = 2.0,
  overbook_t3 = 2.5,
  overbook_t4 = 3.0,

  # --- Deadline expiry behavior ---
  # "accept": rational author accepts best held offer when deadline expires.
  # "lapse":  offer withdrawn by journal; slot freed, author loses the offer.
  # Default "accept" is the more realistic behavioral assumption. "lapse"
  # available for robustness analysis.
  deadline_behavior = "accept",

  # --- Aspiration (separate from submission set) ---
  # Scale factor for aspiration rank: aspiration = 1 + (1-prestige)*J*scale + noise.
  # Baseline 0.4 maps mean-prestige (0.286) to aspiration ~58.
  # Sensitivity: 0.2 (optimistic authors) to 0.6 (pessimistic authors).
  aspiration_scale = 0.4,

  # --- Simulation ---
  T_steps     = 60,         # Time steps in simulation (one step = one day;
                            # must exceed Scholastica's 47-day avg rejection)

  # --- Monte Carlo ---
  n_sims      = 2000,       # Number of Monte Carlo replications

  # --- Market tightness ---
  # Multiplier on all slot ranges. 1.0 = baseline calibration (~1545 total slots
  # for N=800 articles, ratio ~1.9:1). Lower values create tighter markets where
  # mechanisms differ more. E.g., 0.5 => ~773 slots (tight), 0.25 => ~386 (very tight).
  slot_scale  = 1.0
)

# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

#' Compute actual Pearson correlation between prestige and quality
#' for a given mixing weight rho.
#'
#' General formula (derivation from prestige = rho*q + (1-rho)*eps):
#'   cor = rho * sd_q / sqrt(rho^2 * var_q + (1-rho)^2 * var_eps)
#'
#' When q ~ Beta(q_alpha, q_beta) and eps ~ Beta(2, 5), this reduces to
#' rho / sqrt(rho^2 + (1-rho)^2) only when Var(q) == Var(eps). The general
#' form handles arbitrary Beta parameters for the quality distribution.
#'
#' @param rho     Mixing weight (scalar or vector)
#' @param q_alpha Shape1 of quality Beta distribution (default 2)
#' @param q_beta  Shape2 of quality Beta distribution (default 5)
#' @param eps_alpha Shape1 of epsilon Beta distribution (default 2)
#' @param eps_beta  Shape2 of epsilon Beta distribution (default 5)
compute_actual_rho <- function(rho, q_alpha = 2, q_beta = 5,
                                eps_alpha = 2, eps_beta = 5) {
  var_q   <- (q_alpha * q_beta) / ((q_alpha + q_beta)^2 * (q_alpha + q_beta + 1))
  var_eps <- (eps_alpha * eps_beta) / ((eps_alpha + eps_beta)^2 * (eps_alpha + eps_beta + 1))
  rho * sqrt(var_q) / sqrt(rho^2 * var_q + (1 - rho)^2 * var_eps)
}

#' Quality-dependent review: determine how many articles a journal can review
#' given its capacity (in effort units). Low-signal articles are rejected
#' quickly (cost ~0.3 effort units), while high-signal articles require full
#' review (cost 1.0 effort units). This captures the empirical finding that
#' editors spend only minutes rejecting poor submissions but much longer on
#' promising ones (Christensen & Oseid, 2007, Figure 3, p. 198).
#'
#' @param queue_sorted  Integer vector of article IDs in review order
#' @param capacity      Total review effort units available this time step
#' @param signal_mat    J x N signal matrix
#' @param jj            Journal index (row in signal_mat)
#' @return Integer: number of articles from queue_sorted that can be reviewed
compute_review_count <- function(queue_sorted, capacity, signal_mat, jj) {
  if (length(queue_sorted) == 0) return(0L)
  if (capacity <= 0) return(0L)
  sigs <- signal_mat[jj, queue_sorted]
  # Effort cost per article: low-signal articles cost less to review.
  # Scale: articles below median signal cost 0.3 effort, above cost 1.0,
  # with a smooth transition. Use a simple logistic transform.
  sig_median <- median(sigs)
  sig_range <- max(sigs) - min(sigs)
  if (sig_range < 1e-10) {
    # All signals roughly equal — each costs 0.5 effort
    effort <- rep(0.5, length(sigs))
  } else {
    # Standardize signal, then apply logistic: ranges from ~0.3 to ~1.0
    z <- (sigs - sig_median) / (sig_range / 4)
    effort <- 0.3 + 0.7 / (1 + exp(-z))
  }
  # Walk through queue, accumulating effort until capacity exhausted
  cum_effort <- cumsum(effort)
  n_review <- sum(cum_effort <= capacity)
  # Always review at least 1 if queue is non-empty and capacity > 0
  max(1L, as.integer(n_review))
}

#' Build submission queues: vectorized construction of (journal, article) pairs.
#' Returns a list of length J, where element jj is an integer vector of
#' article IDs submitted to journal jj.
build_queues <- function(articles, J) {
  lens <- articles$sub_hi - articles$sub_lo + 1L
  j_col <- unlist(Map(seq.int, articles$sub_lo, articles$sub_hi))
  a_col <- rep(seq_len(nrow(articles)), lens)
  queue_map <- split(a_col, j_col)

  full_queue <- vector("list", J)
  for (jj in seq_len(J)) {
    key <- as.character(jj)
    if (key %in% names(queue_map)) {
      q <- queue_map[[key]]
      # Randomize initial queue order within each journal. All articles
      # arrive simultaneously at t=1 (batch submission), so there is no
      # natural arrival order. We break ties by random permutation.
      full_queue[[jj]] <- q[sample.int(length(q))]
    } else {
      full_queue[[jj]] <- integer(0)
    }
  }
  full_queue
}

#' Precompute all random draws needed within a replication for full CRN.
#' This ensures that every mechanism sees the same:
#'   - initial queue orders (from build_queues)
#'   - per-period journal processing order
#'   - per-period, per-journal expedite coin flips
#' Returns a list with components: j_queue, j_orders (T x J matrix),
#' expedite_coins (T x J matrix of U[0,1] draws).
make_crn_draws <- function(articles, p) {
  J <- p$J
  T_steps <- p$T_steps

  # Shared initial queue order
  j_queue <- build_queues(articles, J)

  # Per-period journal processing order (each row is a permutation of 1:J)
  j_orders <- matrix(0L, nrow = T_steps, ncol = J)
  for (t in 1:T_steps) {
    j_orders[t, ] <- sample.int(J)
  }

  # Per-period, per-journal expedite honor coin flips
  expedite_coins <- matrix(runif(T_steps * J), nrow = T_steps, ncol = J)

  list(j_queue = j_queue, j_orders = j_orders, expedite_coins = expedite_coins)
}

#' Precompute J x N signal matrix.
#' signal_mat[j, i] = w_j * quality_i + (1-w_j) * prestige_i + noise_{j,i}
#' Each journal has its own quality-vs-prestige weight w_j (heterogeneity)
#' and independent noise draw.
make_signal_matrix <- function(articles, p, journals = NULL) {
  J <- p$J
  N <- nrow(articles)

  # Use journal-specific w_j if available, otherwise fall back to global w
  if (!is.null(journals) && "w_j" %in% names(journals)) {
    w_vec <- journals$w_j  # length J
  } else {
    w_vec <- rep(p$w, J)
  }

  # Build J x N deterministic component: each row j uses its own w_j
  # w_vec[j] * quality[i] + (1 - w_vec[j]) * prestige[i]
  determ <- outer(w_vec, articles$quality) +
            outer(1 - w_vec, articles$prestige)

  # Add independent noise
  determ + matrix(rnorm(J * N, 0, p$sigma), nrow = J, ncol = N)
}

#' Resolve overbooking and reallocate bumped articles.
#'
#' After the main simulation loop, some journals may hold more articles than
#' their slot count allows. This function: (1) trims each overbooked journal
#' to its slot limit, keeping the highest-signal articles; (2) reallocates
#' bumped articles via a greedy pass — journals with remaining capacity review
#' bumped articles in their submission pool, ranked by signal, and fill slots.
#'
#' @param a_matched_j  Integer vector of length N (journal assignments, NA = unmatched)
#' @param journals     data.table of journal attributes
#' @param articles     data.table of article attributes (needs sub_lo, sub_hi)
#' @param signal_mat   J x N signal matrix
#' @param J            Number of journals
#' @return Updated a_matched_j vector with feasible assignments
resolve_overbooking <- function(a_matched_j, journals, articles, signal_mat, J) {

  # --- Phase 1: Trim overbooked journals ---
  # For each journal with more matched articles than slots, keep the best
  # by signal and bump the rest. Track which articles are bumped so Phase 2
  # only reallocates those specific articles.
  j_capacity_left <- journals$slots  # will track remaining capacity
  bumped_ids <- integer(0)           # collect IDs of articles bumped in Phase 1

  for (jj in 1:J) {
    matched_here <- which(a_matched_j == jj)
    n_here <- length(matched_here)
    if (n_here > journals$slots[jj]) {
      sigs <- signal_mat[jj, matched_here]
      keep_order <- order(sigs, decreasing = TRUE)
      n_keep <- journals$slots[jj]
      drop <- matched_here[keep_order[(n_keep + 1):n_here]]
      a_matched_j[drop] <- NA_integer_
      bumped_ids <- c(bumped_ids, drop)
      j_capacity_left[jj] <- 0L
    } else {
      j_capacity_left[jj] <- journals$slots[jj] - n_here
    }
  }

  # --- Phase 2: Reallocate ONLY articles bumped in Phase 1 ---
  # Greedy pass: process journals in prestige order (rank 1 first).
  # Each journal with remaining capacity reviews bumped articles that
  # submitted to it, ranks by signal, and fills its open slots.
  # IMPORTANT: We only consider articles that were actually bumped above,
  # not all unmatched articles. Articles that were never matched in the
  # main simulation (e.g., never received an offer) remain unmatched.
  if (length(bumped_ids) == 0) return(a_matched_j)

  # Track which bumped articles are still unplaced
  is_bumped <- rep(FALSE, length(a_matched_j))
  is_bumped[bumped_ids] <- TRUE

  for (jj in 1:J) {   # journals in rank order (rank 1 = most prestigious)
    if (j_capacity_left[jj] <= 0L) next

    # Find bumped articles that submitted to this journal
    candidates <- which(is_bumped &
                        articles$sub_lo <= jj &
                        articles$sub_hi >= jj)
    if (length(candidates) == 0) next

    # Rank by signal, take up to remaining capacity
    sigs <- signal_mat[jj, candidates]
    ranked <- candidates[order(sigs, decreasing = TRUE)]
    n_accept <- min(length(ranked), j_capacity_left[jj])
    accepted <- ranked[1:n_accept]

    a_matched_j[accepted] <- jj
    is_bumped[accepted] <- FALSE
    j_capacity_left[jj] <- j_capacity_left[jj] - n_accept
  }

  a_matched_j
}

# ==============================================================================
# 3. AGENT INITIALIZATION
# ==============================================================================

init_journals <- function(p) {
  J <- p$J
  ss <- if (!is.null(p$slot_scale)) p$slot_scale else 1.0

  tier <- ifelse(1:J <= 20, 1,
           ifelse(1:J <= 50, 2,
            ifelse(1:J <= 100, 3, 4)))

  slots <- integer(J)
  capacity <- integer(J)
  deadline <- integer(J)
  overbook <- numeric(J)
  expedite_prob <- numeric(J)

  # Expedite responsiveness: journal-by-journal draw from Beta, no tier gradient
  ep_a <- p$expedite_prob_mean * p$expedite_prob_conc
  ep_b <- (1 - p$expedite_prob_mean) * p$expedite_prob_conc

  for (i in 1:J) {
    # Expedite prob is tier-independent (drawn per journal)
    expedite_prob[i] <- rbeta(1, ep_a, ep_b)

    if (tier[i] == 1) {
      slots[i]    <- max(1L, round(sample(p$slots_t1[1]:p$slots_t1[2], 1) * ss))
      capacity[i] <- sample(p$cap_t1[1]:p$cap_t1[2], 1)
      deadline[i] <- sample(p$deadline_t1, 1)
      overbook[i] <- p$overbook_t1
    } else if (tier[i] == 2) {
      slots[i]    <- max(1L, round(sample(p$slots_t2[1]:p$slots_t2[2], 1) * ss))
      capacity[i] <- sample(p$cap_t2[1]:p$cap_t2[2], 1)
      deadline[i] <- sample(p$deadline_other, 1)
      overbook[i] <- p$overbook_t2
    } else if (tier[i] == 3) {
      slots[i]    <- max(1L, round(sample(p$slots_t3[1]:p$slots_t3[2], 1) * ss))
      capacity[i] <- sample(p$cap_t3[1]:p$cap_t3[2], 1)
      deadline[i] <- sample(p$deadline_other, 1)
      overbook[i] <- p$overbook_t3
    } else {
      slots[i]    <- max(1L, round(sample(p$slots_t4[1]:p$slots_t4[2], 1) * ss))
      capacity[i] <- sample(p$cap_t4[1]:p$cap_t4[2], 1)
      deadline[i] <- sample(p$deadline_other, 1)
      overbook[i] <- p$overbook_t4
    }
  }

  # --- Journal-specific signal weights (heterogeneity in how much each
  #     journal relies on credentials vs. article quality) ---
  # Draw w_j from a Beta distribution centered near the global w, with
  # tier-dependent means: top journals rely MORE on credentials (lower w_j),
  # lower journals rely more on quality signals (higher w_j).
  # This matches Christensen & Oseid (2007, Tables 1-4): 83-100% of Top 25
  # editors are influenced by credentials vs. 44-56% at 4th Tier.
  w_j <- numeric(J)
  for (i in 1:J) {
    # Tier-dependent mean for w_j (quality weight):
    #   Tier 1: lower w (more prestige-reliant) ~0.50
    #   Tier 2: moderate ~0.55
    #   Tier 3: moderate-high ~0.65
    #   Tier 4: higher w (more quality-reliant) ~0.70
    w_mean <- if (tier[i] == 1) 0.50
              else if (tier[i] == 2) 0.55
              else if (tier[i] == 3) 0.65
              else 0.70
    # Beta params with moderate spread (sd ~ 0.08)
    w_conc <- 25  # concentration parameter
    w_a <- w_mean * w_conc
    w_b <- (1 - w_mean) * w_conc
    w_j[i] <- rbeta(1, w_a, w_b)
  }

  data.table(
    j_id     = 1:J,
    rank     = 1:J,
    tier     = tier,
    slots    = slots,
    capacity = capacity,
    deadline = deadline,
    overbook = overbook,
    expedite_prob = expedite_prob,
    w_j      = w_j,
    prestige = 1 - (0:(J-1)) / (J - 1)   # rank 1 -> 1.0, rank J -> 0.0
  )
}

init_articles <- function(p) {
  N <- p$N
  J <- p$J

  quality  <- rbeta(N, p$q_alpha, p$q_beta)
  epsilon  <- rbeta(N, p$eps_alpha, p$eps_beta)
  # prestige is a mixture of quality and independent noise; rho is the mixing
  # weight (see compute_actual_rho for the implied Pearson correlation).
  prestige <- pmin(1, pmax(0, p$rho * quality + (1 - p$rho) * epsilon))

  # --- Top-heavy submission sets ---
  # In reality, the near-zero marginal cost of electronic submission (Scholastica,
  # formerly ExpressO) means almost every author submits to elite journals, and
  # the main variation is how far DOWN the hierarchy an author also submits.
  # Georgetown (rank ~25) receives over 2,000 submissions per year; elite journals
  # are flooded. This produces the empirically realistic pattern of heavy
  # congestion at the top and lighter loads at the bottom.
  #
  # sub_lo: nearly everyone starts near rank 1.  The marginal cost of adding
  # Harvard or Yale to a Scholastica submission batch is essentially zero, so
  # sub_lo is not prestige-dependent — it is centered at 1 with small noise.
  #   ~55% of authors have sub_lo = 1
  #   ~82% have sub_lo <= 5
  #   ~97% have sub_lo <= 10
  sub_lo <- pmax(1L, as.integer(round(rnorm(N, 1, 5))))
  sub_lo <- pmin(sub_lo, as.integer(J))

  # sub_hi: prestige-dependent reach downward, using an EXPONENTIAL
  # distribution for width. The exponential naturally produces a monotonically
  # declining submission density curve: P(width >= r) = exp(-r / mean_width),
  # so higher-ranked (more prestigious) journals get more submissions and
  # lower-ranked journals get progressively fewer. This avoids the plateau
  # artifact that arises from deterministic or uniform-width formulas.
  #
  # Mean width depends on prestige:
  #   High-prestige authors (p~0.8): mean ~40 (targeted, elite-focused)
  #   Mean-prestige authors (p~0.3): mean ~65 (moderate spread)
  #   Low-prestige authors (p~0.1): mean ~75 (wide net)
  mean_width <- 30 + (1 - prestige) * 50
  sub_width <- pmax(10L, as.integer(round(rexp(N, rate = 1 / mean_width))))

  sub_hi <- pmin(as.integer(J), sub_lo + sub_width - 1L)
  # Final safety: ensure sub_lo <= sub_hi
  sub_lo <- pmin(sub_lo, sub_hi)

  # --- Aspiration rank (separate from submission set) ---
  # Where an author submits (sub_lo near 1 for everyone) and where an
  # author realistically aspires to place are different things. Adding
  # Harvard to a Scholastica batch is costless; genuinely expecting
  # Harvard to accept is another matter. We define an aspiration rank
  # a_i that reflects an author's realistic self-assessment, anchored to
  # prestige with noise:
  #   High-prestige (p~0.8): aspiration ~17 (expects top-20 placement)
  #   Mean-prestige (p~0.3): aspiration ~57 (expects mid-market)
  #   Low-prestige (p~0.1): aspiration ~73 (expects lower half)
  # The acceptance threshold (Section 4) uses aspiration_rank, NOT sub_lo,
  # to determine when an author accepts an offer vs. holds and expedites.
  # Compute unconstrained latent aspiration, then clamp to [sub_lo, sub_hi].
  # Without clamping:
  #   - If a_i > h_i, the threshold formula inverts (authors become MORE selective
  #     over time, which reverses the intended behavior).
  #   - If a_i < ℓ_i, the author "aspires" to journals they never submitted to,
  #     which is conceptually incoherent given near-zero marginal submission costs.
  aspiration_latent <- as.integer(round(
    1 + (1 - prestige) * J * p$aspiration_scale + rnorm(N, 0, 10)
  ))
  aspiration_rank <- pmin(pmax(aspiration_latent, sub_lo), sub_hi)

  data.table(
    a_id     = 1:N,
    quality  = quality,
    prestige = prestige,
    sub_lo   = sub_lo,
    sub_hi   = sub_hi,
    aspiration = aspiration_rank
  )
}

# ==============================================================================
# 4. CURRENT SYSTEM SIMULATION (Decentralized Expedite)
# ==============================================================================

simulate_current <- function(journals, articles, p, signal_mat = NULL,
                             crn = NULL) {
  J <- p$J
  N <- p$N
  T_steps <- p$T_steps
  deadline_accept <- identical(p$deadline_behavior, "accept")

  # Use shared signal matrix if provided (common random numbers),
  # otherwise generate an independent one
  if (is.null(signal_mat)) signal_mat <- make_signal_matrix(articles, p, journals)

  # Use shared CRN draws if provided, otherwise generate independently
  if (is.null(crn)) crn <- make_crn_draws(articles, p)
  j_queue <- lapply(crn$j_queue, function(x) x)  # deep copy queues

  # State tracking
  j_remaining     <- journals$slots
  a_best_offer    <- rep(Inf, N)         # best (lowest rank) offer held
  a_best_j        <- rep(NA_integer_, N) # journal id of best offer
  a_matched_j     <- rep(NA_integer_, N) # final matched journal
  a_status        <- rep(STATUS_SUBMITTED, N)
  a_offer_expires <- rep(Inf, N)         # time step when current offer expires

  for (t in 1:T_steps) {

    # --- Handle expired offers ---
    expired <- which(a_offer_expires <= t & a_status == STATUS_HOLDING)
    if (length(expired) > 0) {
      if (deadline_accept) {
        # Rational author accepts their best held offer when deadline forces it
        a_matched_j[expired] <- a_best_j[expired]
        a_status[expired] <- STATUS_MATCHED
      } else {
        # "Lapse" behavior: journal withdraws offer, slot freed
        for (ei in expired) {
          j_remaining[a_best_j[ei]] <- j_remaining[a_best_j[ei]] + 1L
          a_best_offer[ei] <- Inf
          a_best_j[ei] <- NA_integer_
          a_status[ei] <- STATUS_SUBMITTED
          a_offer_expires[ei] <- Inf
        }
      }
    }

    # Use precomputed journal processing order from CRN draws
    j_order <- crn$j_orders[t, ]

    for (jj in j_order) {
      # A journal can issue offers up to its overbooking cap, not just its
      # actual slot count. j_committed = slots - j_remaining counts offers
      # currently outstanding or accepted. The journal stops only when
      # committed >= ceil(slots * overbook).
      j_committed <- journals$slots[jj] - j_remaining[jj]
      overbook_cap <- ceiling(journals$slots[jj] * journals$overbook[jj])
      if (j_committed >= overbook_cap) next
      if (length(j_queue[[jj]]) == 0) next

      # Remove already-matched articles from queue
      queue <- j_queue[[jj]]
      queue <- queue[a_status[queue] != STATUS_MATCHED]
      j_queue[[jj]] <- queue
      if (length(queue) == 0) next

      # Expedite prioritization: articles holding offers from other journals
      # are moved to the front of the queue, BUT only if this journal honors
      # the expedite request (probabilistic, journal-level heterogeneity).
      # This captures the real-world variation: some journals honor expedites
      # routinely, others say expedites confer no competitive advantage.
      is_expedite <- a_status[queue] == STATUS_HOLDING
      if (any(is_expedite) && crn$expedite_coins[t, jj] <= journals$expedite_prob[jj]) {
        queue_sorted <- c(queue[is_expedite], queue[!is_expedite])
      } else {
        queue_sorted <- queue  # no expedite priority this round
      }

      # Review with quality-dependent effort: bad articles are rejected
      # quickly, good articles take more review time
      n_review <- compute_review_count(queue_sorted, journals$capacity[jj],
                                        signal_mat, jj)
      n_review <- min(n_review, length(queue_sorted))
      reviewed <- queue_sorted[1:n_review]

      # Remove reviewed articles from queue (each article evaluated once)
      j_queue[[jj]] <- queue_sorted[-(1:n_review)]

      # Signals from precomputed matrix
      sigs <- signal_mat[jj, reviewed]

      # Number of new offers: overbooking cap minus already-committed slots
      n_offers <- min(overbook_cap - j_committed, length(reviewed))
      n_offers <- max(n_offers, 0L)
      if (n_offers == 0) next

      top_idx <- order(sigs, decreasing = TRUE)[1:n_offers]
      offer_articles <- reviewed[top_idx]

      for (ai in offer_articles) {
        if (a_status[ai] == STATUS_MATCHED) next

        # Is this offer better than the author's current best?
        # (When a_best_j is NA, a_best_offer is Inf, so any finite jj wins.)
        if (jj < a_best_offer[ai]) {
          old_j <- a_best_j[ai]

          a_best_offer[ai] <- jj
          a_best_j[ai] <- jj
          a_status[ai] <- STATUS_HOLDING
          a_offer_expires[ai] <- t + journals$deadline[jj]

          # Return slot to the previous journal (trading up)
          if (!is.na(old_j)) {
            j_remaining[old_j] <- j_remaining[old_j] + 1L
          }
          # Reserve a slot at the new journal
          j_remaining[jj] <- j_remaining[jj] - 1L

          # --- Time-dependent acceptance threshold ---
          # Uses aspiration rank (not sub_lo) as the ceiling. Submitting
          # to Harvard is costless; genuinely aspiring to Harvard is not.
          # Early: accept if within top 10% of aspiration-to-sub_hi range.
          # Late: accept top 60%. Captures cascade strategy (Dorf 2012).
          time_frac <- t / T_steps  # 0 at start, 1 at end
          # Acceptance threshold rises over time (more of the range is acceptable)
          accept_pct <- 0.10 + 0.50 * time_frac  # 10% early -> 60% late
          top_threshold <- articles$aspiration[ai] +
            accept_pct * (articles$sub_hi[ai] - articles$aspiration[ai])
          if (jj <= floor(top_threshold)) {
            a_matched_j[ai] <- jj
            a_status[ai] <- STATUS_MATCHED
          }
        }
        # If offer is worse than current best, author ignores it
      }
    }
  }

  # Finalize: anyone still holding accepts their best offer
  still_holding <- which(a_status == STATUS_HOLDING)
  a_matched_j[still_holding] <- a_best_j[still_holding]

  # --- Overbooking resolution + reallocation ---
  a_matched_j <- resolve_overbooking(a_matched_j, journals, articles, signal_mat, J)

  data.table(
    a_id      = 1:N,
    quality   = articles$quality,
    prestige  = articles$prestige,
    matched_j = a_matched_j
  )
}

# ==============================================================================
# 5. DEFERRED ACCEPTANCE (Centralized Match)
# ==============================================================================
# Each journal forms its own noisy preference list (heterogeneous signals),
# matching the information structure of the decentralized mechanisms.
# This ensures the DA vs. current-system comparison isolates mechanism design,
# not information advantage.

simulate_da <- function(journals, articles, p, proposer = "journal",
                         signal_mat = NULL) {
  J <- p$J
  N <- p$N

  # Use shared signal matrix if provided (common random numbers),
  # otherwise generate an independent one
  if (is.null(signal_mat)) signal_mat <- make_signal_matrix(articles, p, journals)

  if (proposer == "journal") {
    # Journal preference lists: each journal ranks the articles that submitted
    # to it, ordered by that journal's own noisy signal (descending)
    j_prefs <- lapply(1:J, function(jj) {
      submitted <- which(articles$sub_lo <= jj & articles$sub_hi >= jj)
      if (length(submitted) == 0) return(integer(0))
      submitted[order(signal_mat[jj, submitted], decreasing = TRUE)]
    })

    # Cache preference list lengths
    j_pref_lengths <- vapply(j_prefs, length, integer(1))

    # State
    j_slots_left <- journals$slots
    j_next_proposal <- rep(1L, J)
    a_held_by <- rep(NA_integer_, N)
    a_held_rank <- rep(Inf, N)

    max_rounds <- N * J
    round <- 0

    # Batch-propose: each journal proposes to enough candidates per round
    # to potentially fill all its open slots. This is an optimization that
    # converges to the same stable matching as one-at-a-time DA.
    while (round < max_rounds) {
      round <- round + 1
      made_proposal <- FALSE

      free_journals <- which(j_slots_left > 0 & j_next_proposal <= j_pref_lengths)
      if (length(free_journals) == 0) break

      for (jj in free_journals) {
        if (j_slots_left[jj] <= 0L) next
        plist <- j_prefs[[jj]]
        start_idx <- j_next_proposal[jj]
        if (start_idx > j_pref_lengths[jj]) next

        end_idx <- min(start_idx + j_slots_left[jj] - 1L, j_pref_lengths[jj])
        candidates <- plist[start_idx:end_idx]
        j_next_proposal[jj] <- end_idx + 1L
        made_proposal <- TRUE

        for (ai in candidates) {
          # Article prefers whichever journal has lower rank (= more prestigious)
          if (jj < a_held_rank[ai]) {
            old_j <- a_held_by[ai]
            if (!is.na(old_j)) {
              j_slots_left[old_j] <- j_slots_left[old_j] + 1L
            }
            a_held_by[ai] <- jj
            a_held_rank[ai] <- jj
            j_slots_left[jj] <- j_slots_left[jj] - 1L
          }
        }
      }

      if (!made_proposal) break
    }

    return(data.table(
      a_id      = 1:N,
      quality   = articles$quality,
      prestige  = articles$prestige,
      matched_j = a_held_by
    ))

  } else if (proposer == "journal_constrained") {
    # ---- Capacity-constrained journal-proposing DA ----
    # Each journal can only rank the top K articles it could feasibly review
    # given its per-step capacity over T time steps. This prevents the
    # unconstrained DA from giving journals an unrealistic information
    # advantage (ranking articles they'd never have time to read).
    #
    # K_j = floor(capacity_j * T_steps / mean_effort), where mean_effort
    # is computed from the actual signal distribution for that journal's
    # submission pool (matching compute_review_count's effort model).
    j_prefs <- lapply(1:J, function(jj) {
      submitted <- which(articles$sub_lo <= jj & articles$sub_hi >= jj)
      if (length(submitted) == 0) return(integer(0))
      # Rank by signal (same as unconstrained DA)
      ranked <- submitted[order(signal_mat[jj, submitted], decreasing = TRUE)]
      # Compute capacity constraint: how many articles could this journal
      # review over the full T-step simulation?
      total_capacity <- journals$capacity[jj] * p$T_steps
      # Estimate mean effort using compute_review_count's logistic model
      # applied to the ranked articles. Walk through the list accumulating
      # effort until the total capacity is exhausted.
      sigs <- signal_mat[jj, ranked]
      sig_median <- median(sigs)
      sig_range <- max(sigs) - min(sigs)
      if (sig_range < 1e-10) {
        effort <- rep(0.5, length(sigs))
      } else {
        z <- (sigs - sig_median) / (sig_range / 4)
        effort <- 0.3 + 0.7 / (1 + exp(-z))
      }
      cum_effort <- cumsum(effort)
      K <- max(1L, sum(cum_effort <= total_capacity))
      # Truncate preference list to top-K reviewable articles
      ranked[1:min(K, length(ranked))]
    })

    # Cache preference list lengths
    j_pref_lengths <- vapply(j_prefs, length, integer(1))

    # State
    j_slots_left <- journals$slots
    j_next_proposal <- rep(1L, J)
    a_held_by <- rep(NA_integer_, N)
    a_held_rank <- rep(Inf, N)

    max_rounds <- N * J
    round <- 0

    while (round < max_rounds) {
      round <- round + 1
      made_proposal <- FALSE

      free_journals <- which(j_slots_left > 0 & j_next_proposal <= j_pref_lengths)
      if (length(free_journals) == 0) break

      for (jj in free_journals) {
        if (j_slots_left[jj] <= 0L) next
        plist <- j_prefs[[jj]]
        start_idx <- j_next_proposal[jj]
        if (start_idx > j_pref_lengths[jj]) next

        end_idx <- min(start_idx + j_slots_left[jj] - 1L, j_pref_lengths[jj])
        candidates <- plist[start_idx:end_idx]
        j_next_proposal[jj] <- end_idx + 1L
        made_proposal <- TRUE

        for (ai in candidates) {
          if (jj < a_held_rank[ai]) {
            old_j <- a_held_by[ai]
            if (!is.na(old_j)) {
              j_slots_left[old_j] <- j_slots_left[old_j] + 1L
            }
            a_held_by[ai] <- jj
            a_held_rank[ai] <- jj
            j_slots_left[jj] <- j_slots_left[jj] - 1L
          }
        }
      }

      if (!made_proposal) break
    }

    return(data.table(
      a_id      = 1:N,
      quality   = articles$quality,
      prestige  = articles$prestige,
      matched_j = a_held_by
    ))

  } else if (proposer == "author") {
    # ---- Author-proposing DA: author-optimal, journal-pessimal ----
    # Each author proposes to journals in order of prestige (most prestigious
    # first). Journals tentatively hold the best `slots` candidates by signal
    # and reject the rest. The difference between journal- and author-optimal
    # stable matchings bounds the lattice of stable matchings; if they are
    # nearly identical, the stable matching is effectively unique.

    # Author preference lists: journals in submission range, ascending rank
    a_prefs <- lapply(1:N, function(ai) {
      seq.int(articles$sub_lo[ai], articles$sub_hi[ai])
    })
    a_pref_lengths <- vapply(a_prefs, length, integer(1))
    a_next_proposal <- rep(1L, N)

    # Journal holding state
    j_held <- vector("list", J)
    for (jj in 1:J) j_held[[jj]] <- integer(0)

    # Track which journal currently holds each author (NA if free)
    a_held_by <- rep(NA_integer_, N)

    max_rounds <- sum(a_pref_lengths)
    round <- 0

    while (round < max_rounds) {
      round <- round + 1

      # Free authors: not held and still have proposals to make
      free_authors <- which(is.na(a_held_by) &
                            a_next_proposal <= a_pref_lengths)
      if (length(free_authors) == 0) break

      # Each free author proposes to their next-choice journal
      prop_journals <- vapply(free_authors, function(ai) {
        a_prefs[[ai]][a_next_proposal[ai]]
      }, integer(1))
      a_next_proposal[free_authors] <- a_next_proposal[free_authors] + 1L

      # Process proposals by journal: combine held + new, keep best `slots`
      for (target_j in unique(prop_journals)) {
        new_applicants <- free_authors[prop_journals == target_j]
        all_candidates <- c(j_held[[target_j]], new_applicants)

        if (length(all_candidates) <= journals$slots[target_j]) {
          # Room for everyone — accept all
          j_held[[target_j]] <- all_candidates
          a_held_by[all_candidates] <- target_j
        } else {
          # Keep the best `slots` by signal; reject the rest
          sigs <- signal_mat[target_j, all_candidates]
          keep_order <- order(sigs, decreasing = TRUE)
          n_keep <- journals$slots[target_j]
          kept     <- all_candidates[keep_order[1:n_keep]]
          rejected <- all_candidates[keep_order[(n_keep + 1):length(all_candidates)]]

          j_held[[target_j]] <- kept
          a_held_by[kept] <- target_j
          a_held_by[rejected] <- NA_integer_
        }
      }
    }

    return(data.table(
      a_id      = 1:N,
      quality   = articles$quality,
      prestige  = articles$prestige,
      matched_j = a_held_by
    ))

  } else {
    stop("proposer must be 'journal' or 'author', got: ", proposer)
  }
}

# ==============================================================================
# 6. COUNTERFACTUAL: EXTENDED DEADLINES
# ==============================================================================

simulate_extended_deadlines <- function(journals, articles, p,
                                         uniform_deadline = 14,
                                         signal_mat = NULL,
                                         crn = NULL) {
  j_mod <- copy(journals)
  j_mod$deadline <- uniform_deadline
  simulate_current(j_mod, articles, p, signal_mat = signal_mat, crn = crn)
}

# ==============================================================================
# 7. COUNTERFACTUAL: NO EXPEDITE
# ==============================================================================
# Identical to the current system EXCEPT that articles holding offers from
# other journals receive NO queue priority. Journals review in submission
# order regardless of expedite status. Authors can still hold offers and
# trade up — the only change is that holding an offer does not move them
# to the front of other journals' queues.

simulate_no_expedite <- function(journals, articles, p, signal_mat = NULL,
                                  crn = NULL) {
  J <- p$J
  N <- p$N
  T_steps <- p$T_steps
  deadline_accept <- identical(p$deadline_behavior, "accept")

  if (is.null(signal_mat)) signal_mat <- make_signal_matrix(articles, p, journals)
  if (is.null(crn)) crn <- make_crn_draws(articles, p)
  j_queue <- lapply(crn$j_queue, function(x) x)  # deep copy queues

  j_remaining     <- journals$slots
  a_best_offer    <- rep(Inf, N)
  a_best_j        <- rep(NA_integer_, N)
  a_matched_j     <- rep(NA_integer_, N)
  a_status        <- rep(STATUS_SUBMITTED, N)
  a_offer_expires <- rep(Inf, N)

  for (t in 1:T_steps) {

    # Handle expired offers (same logic as current system)
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
      # Overbooking cap check (same logic as simulate_current)
      j_committed <- journals$slots[jj] - j_remaining[jj]
      overbook_cap <- ceiling(journals$slots[jj] * journals$overbook[jj])
      if (j_committed >= overbook_cap) next
      queue <- j_queue[[jj]]
      queue <- queue[a_status[queue] != STATUS_MATCHED]
      j_queue[[jj]] <- queue
      if (length(queue) == 0) next

      # KEY DIFFERENCE: No expedite priority. Review in original queue order
      # (no sorting by holding status). This is the ONLY change from the
      # current system — holding and trading up are still permitted.

      n_review <- compute_review_count(queue, journals$capacity[jj],
                                        signal_mat, jj)
      n_review <- min(n_review, length(queue))
      reviewed <- queue[1:n_review]
      j_queue[[jj]] <- queue[-(1:n_review)]

      sigs <- signal_mat[jj, reviewed]

      n_offers <- min(overbook_cap - j_committed, length(reviewed))
      n_offers <- max(n_offers, 0L)
      if (n_offers == 0) next

      top_idx <- order(sigs, decreasing = TRUE)[1:n_offers]
      offer_articles <- reviewed[top_idx]

      for (ai in offer_articles) {
        if (a_status[ai] == STATUS_MATCHED) next

        # Same offer logic as current system — authors can hold and trade up
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

          # Time-dependent acceptance threshold (same as current system)
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

  # --- Overbooking resolution + reallocation ---
  a_matched_j <- resolve_overbooking(a_matched_j, journals, articles, signal_mat, J)

  data.table(
    a_id      = 1:N,
    quality   = articles$quality,
    prestige  = articles$prestige,
    matched_j = a_matched_j
  )
}

# ==============================================================================
# 7b. COUNTERFACTUAL: NO EXPEDITE + EXTENDED DEADLINES
# ==============================================================================
# Combines no-expedite queue-jumping with uniform 14-day deadlines.
# This isolates decentralization itself: no queue-jumping, no deadline
# pressure, but still a decentralized sequential process. Used in the
# efficiency loss decomposition to ensure each component changes exactly
# one feature at a time.

simulate_no_expedite_extended <- function(journals, articles, p,
                                          uniform_deadline = 14,
                                          signal_mat = NULL,
                                          crn = NULL) {
  j_mod <- copy(journals)
  j_mod$deadline <- uniform_deadline
  simulate_no_expedite(j_mod, articles, p, signal_mat = signal_mat, crn = crn)
}

# ==============================================================================
# 8. COUNTERFACTUAL: NO TRADING UP
# ==============================================================================
# Authors accept the first offer they receive. No holding, no trading up.
# This is a more aggressive restriction than no-expedite and provides a
# distinct counterfactual: what happens when authors cannot strategically
# delay acceptance?

simulate_no_trading_up <- function(journals, articles, p, signal_mat = NULL,
                                    crn = NULL) {
  J <- p$J
  N <- p$N
  T_steps <- p$T_steps

  if (is.null(signal_mat)) signal_mat <- make_signal_matrix(articles, p, journals)
  if (is.null(crn)) crn <- make_crn_draws(articles, p)
  j_queue <- lapply(crn$j_queue, function(x) x)  # deep copy queues

  j_remaining <- journals$slots
  a_matched_j <- rep(NA_integer_, N)
  a_status    <- rep(STATUS_SUBMITTED, N)

  for (t in 1:T_steps) {
    j_order <- crn$j_orders[t, ]

    for (jj in j_order) {
      # Overbooking cap check (same logic as other mechanisms)
      j_committed <- journals$slots[jj] - j_remaining[jj]
      overbook_cap <- ceiling(journals$slots[jj] * journals$overbook[jj])
      if (j_committed >= overbook_cap) next
      queue <- j_queue[[jj]]
      queue <- queue[a_status[queue] != STATUS_MATCHED]
      j_queue[[jj]] <- queue
      if (length(queue) == 0) next

      # Review in submission order (no expedite priority possible)
      n_review <- compute_review_count(queue, journals$capacity[jj],
                                        signal_mat, jj)
      n_review <- min(n_review, length(queue))
      reviewed <- queue[1:n_review]
      j_queue[[jj]] <- queue[-(1:n_review)]

      sigs <- signal_mat[jj, reviewed]

      n_offers <- min(overbook_cap - j_committed, length(reviewed))
      n_offers <- max(n_offers, 0L)
      if (n_offers == 0) next

      top_idx <- order(sigs, decreasing = TRUE)[1:n_offers]
      offer_articles <- reviewed[top_idx]

      for (ai in offer_articles) {
        if (a_status[ai] == STATUS_MATCHED) next

        # Accept first offer unconditionally (no holding, no trading up).
        # Journal-side behavior is unchanged: the journal issues the full
        # batch of overbooked offers this period, just as in the current system.
        a_matched_j[ai] <- jj
        a_status[ai] <- STATUS_MATCHED
        j_remaining[jj] <- j_remaining[jj] - 1L
      }
    }
  }

  # --- Overbooking resolution + reallocation ---
  # Same as other sequential mechanisms: journals may hold more articles
  # than slots due to overbooking, so we trim and reallocate.
  a_matched_j <- resolve_overbooking(a_matched_j, journals, articles, signal_mat, J)

  data.table(
    a_id      = 1:N,
    quality   = articles$quality,
    prestige  = articles$prestige,
    matched_j = a_matched_j
  )
}

# ==============================================================================
# 9. COUNTERFACTUAL: RANDOM ASSIGNMENT (lower bound)
# ==============================================================================

simulate_random <- function(journals, articles, p) {
  N <- p$N
  slot_pool <- rep(journals$j_id, journals$slots)

  if (length(slot_pool) >= N) {
    assigned_slots <- sample(slot_pool, N, replace = FALSE)
  } else {
    assigned_slots <- rep(NA_integer_, N)
    lucky <- sample(N, length(slot_pool))
    assigned_slots[lucky] <- sample(slot_pool)
  }

  data.table(
    a_id      = 1:N,
    quality   = articles$quality,
    prestige  = articles$prestige,
    matched_j = assigned_slots
  )
}

# ==============================================================================
# 10. WELFARE MEASURES
# ==============================================================================
# Takes the full articles table (for stable quartile computation) plus the
# matching result and journal data.

compute_welfare <- function(result, journals, articles_full, p) {
  J <- p$J
  N <- p$N

  matched <- copy(result[!is.na(matched_j)])
  n_matched <- nrow(matched)
  n_unmatched <- N - n_matched

  if (n_matched == 0) {
    return(list(
      match_quality = 0,
      match_quality_allN = 0,
      rank_cor = 0,
      n_matched = 0,
      n_unmatched = N,
      tier_quality = rep(0, 4),
      tier_avg_article_q = rep(0, 4),
      tier_var_article_q = rep(0, 4),
      tier_n = rep(0L, 4),
      author_avg_q_by_pquartile = rep(0, 4),
      author_avg_jrank_by_pquartile = rep(0, 4),
      cond_prestige_bias = lapply(1:5, function(x)
        list(mean_jrank = rep(NA_real_, 4), n = rep(0L, 4)))
    ))
  }

  # Journal prestige score
  j_prestige <- 1 - (matched$matched_j - 1) / (J - 1)

  # Aggregate match quality (supermodular payoff)
  agg_quality <- mean(matched$quality * j_prestige)

  # All-N match quality: total payoff / N (unmatched articles contribute 0).
  # This penalizes mechanisms that leave articles unmatched, enabling fair
  # comparison across mechanisms with different match rates.
  agg_quality_allN <- sum(matched$quality * j_prestige) / N

  # Rank correlation
  rank_cor <- cor(matched$quality, j_prestige, method = "spearman")

  # --- Tier-level analysis ---
  matched[, j_tier := ifelse(matched_j <= 20, 1,
                       ifelse(matched_j <= 50, 2,
                        ifelse(matched_j <= 100, 3, 4)))]

  tier_quality <- numeric(4)
  tier_avg_q   <- numeric(4)
  tier_var_q   <- numeric(4)
  tier_n       <- integer(4)

  for (tt in 1:4) {
    tier_rows <- matched[j_tier == tt]
    if (nrow(tier_rows) > 0) {
      jp <- 1 - (tier_rows$matched_j - 1) / (J - 1)
      tier_quality[tt] <- mean(tier_rows$quality * jp)
      tier_avg_q[tt]   <- mean(tier_rows$quality)
      tier_var_q[tt]   <- if (nrow(tier_rows) >= 2) var(tier_rows$quality) else 0
      tier_n[tt]       <- nrow(tier_rows)
    }
  }

  # --- Author prestige quartile analysis ---
  # Quartile breaks are computed on the FULL population so that bins are
  # stable across mechanisms (even if different mechanisms leave different
  # articles unmatched).
  full_breaks <- quantile(articles_full$prestige, probs = 0:4/4, na.rm = TRUE)
  matched[, a_quartile := cut(prestige,
    breaks = full_breaks, labels = 1:4, include.lowest = TRUE)]

  author_q_by_quartile <- numeric(4)
  author_j_rank_by_quartile <- numeric(4)

  for (qq in 1:4) {
    q_rows <- matched[a_quartile == qq]
    if (nrow(q_rows) > 0) {
      author_q_by_quartile[qq] <- mean(q_rows$quality)
      author_j_rank_by_quartile[qq] <- mean(q_rows$matched_j)
    }
  }

  # --- Conditional prestige-bias analysis ---
  # Bin articles by TRUE quality quintile (computed on full population),
  # then within each quintile report mean placement rank by prestige quartile.
  # This isolates prestige bias: if two articles of equal quality get different
  # placements based on prestige, that is inequity.
  full_q_breaks <- quantile(articles_full$quality, probs = 0:5/5, na.rm = TRUE)
  # Ensure strictly increasing breaks
  full_q_breaks <- unique(full_q_breaks)
  n_qbins <- length(full_q_breaks) - 1

  # Build quality quintile labels for the full population
  result_copy <- copy(result)
  result_copy[, q_quintile := cut(quality,
    breaks = full_q_breaks, labels = seq_len(n_qbins), include.lowest = TRUE)]
  result_copy[, p_quartile := cut(prestige,
    breaks = quantile(articles_full$prestige, probs = 0:4/4, na.rm = TRUE),
    labels = 1:4, include.lowest = TRUE)]

  # For matched articles, compute mean journal rank by (quality quintile, prestige quartile)
  matched_copy <- result_copy[!is.na(matched_j)]
  cond_prestige_bias <- list()
  for (qq in seq_len(n_qbins)) {
    row <- numeric(4)
    n_row <- integer(4)
    for (pp in 1:4) {
      sub <- matched_copy[q_quintile == qq & p_quartile == pp]
      if (nrow(sub) > 0) {
        row[pp] <- mean(sub$matched_j)
        n_row[pp] <- nrow(sub)
      } else {
        row[pp] <- NA_real_
        n_row[pp] <- 0L
      }
    }
    cond_prestige_bias[[qq]] <- list(mean_jrank = row, n = n_row)
  }

  list(
    match_quality         = agg_quality,
    match_quality_allN    = agg_quality_allN,
    rank_cor              = rank_cor,
    n_matched             = n_matched,
    n_unmatched           = n_unmatched,
    tier_quality          = tier_quality,
    tier_avg_article_q    = tier_avg_q,
    tier_var_article_q    = tier_var_q,
    tier_n                = tier_n,
    author_avg_q_by_pquartile = author_q_by_quartile,
    author_avg_jrank_by_pquartile = author_j_rank_by_quartile,
    cond_prestige_bias = cond_prestige_bias
  )
}

# ==============================================================================
# 11. SINGLE REPLICATION (runs all mechanisms on same population)
# ==============================================================================

run_one_replication <- function(sim_id, p) {
  journals <- init_journals(p)
  articles <- init_articles(p)

  # Common random numbers: all mechanisms see the same editorial signals,
  # initial queue orders, per-period journal processing orders, and
  # expedite coin flips. This ensures mechanism comparisons reflect only
  # design differences, not Monte Carlo noise.
  signal_mat <- make_signal_matrix(articles, p, journals)
  crn <- make_crn_draws(articles, p)

  res_current   <- simulate_current(journals, articles, p, signal_mat = signal_mat,
                                     crn = crn)
  res_da        <- simulate_da(journals, articles, p, proposer = "journal",
                                signal_mat = signal_mat)
  res_da_author <- simulate_da(journals, articles, p, proposer = "author",
                                signal_mat = signal_mat)
  res_extended  <- simulate_extended_deadlines(journals, articles, p,
                                                signal_mat = signal_mat,
                                                crn = crn)
  res_noexp     <- simulate_no_expedite(journals, articles, p, signal_mat = signal_mat,
                                         crn = crn)
  res_noexp_ext <- simulate_no_expedite_extended(journals, articles, p,
                                                  signal_mat = signal_mat,
                                                  crn = crn)
  res_notradeup <- simulate_no_trading_up(journals, articles, p, signal_mat = signal_mat,
                                           crn = crn)
  res_random    <- simulate_random(journals, articles, p)
  res_da_constr <- simulate_da(journals, articles, p,
                                proposer = "journal_constrained",
                                signal_mat = signal_mat)

  # Pass full articles for stable quartile computation
  w_current    <- compute_welfare(res_current,   journals, articles, p)
  w_da         <- compute_welfare(res_da,        journals, articles, p)
  w_da_author  <- compute_welfare(res_da_author, journals, articles, p)
  w_da_constr  <- compute_welfare(res_da_constr, journals, articles, p)
  w_extended   <- compute_welfare(res_extended,  journals, articles, p)
  w_noexp      <- compute_welfare(res_noexp,     journals, articles, p)
  w_noexp_ext  <- compute_welfare(res_noexp_ext, journals, articles, p)
  w_notradeup  <- compute_welfare(res_notradeup, journals, articles, p)
  w_random     <- compute_welfare(res_random,    journals, articles, p)

  # --- Cross-mechanism comparison: displacement rates vs. journal-DA ---
  # For each mechanism, compute the fraction of articles in each tier that
  # would be placed in a different tier under the DA benchmark.
  get_tier <- function(j_vec) {
    ifelse(is.na(j_vec), NA_integer_,
      ifelse(j_vec <= 20L, 1L,
        ifelse(j_vec <= 50L, 2L,
          ifelse(j_vec <= 100L, 3L, 4L))))
  }

  da_tiers <- get_tier(res_da$matched_j)

  compute_displacement <- function(mech_result) {
    mech_tiers <- get_tier(mech_result$matched_j)
    disp <- numeric(4)
    for (tt in 1:4) {
      in_tier <- which(mech_tiers == tt)
      if (length(in_tier) > 0) {
        # Fraction whose DA tier differs from their mechanism tier
        # NA in da_tiers means unmatched under DA => count as displaced
        disp[tt] <- mean(is.na(da_tiers[in_tier]) | da_tiers[in_tier] != tt)
      }
    }
    disp
  }

  w_current$displacement_vs_da    <- compute_displacement(res_current)
  w_da$displacement_vs_da         <- rep(0, 4)  # DA vs. itself = 0
  w_da_author$displacement_vs_da  <- compute_displacement(res_da_author)
  w_da_constr$displacement_vs_da  <- compute_displacement(res_da_constr)
  w_extended$displacement_vs_da   <- compute_displacement(res_extended)
  w_noexp$displacement_vs_da      <- compute_displacement(res_noexp)
  w_noexp_ext$displacement_vs_da  <- compute_displacement(res_noexp_ext)
  w_notradeup$displacement_vs_da  <- compute_displacement(res_notradeup)
  w_random$displacement_vs_da     <- compute_displacement(res_random)

  list(
    sim_id     = sim_id,
    current    = w_current,
    da         = w_da,
    da_author  = w_da_author,
    da_constr  = w_da_constr,
    extended   = w_extended,
    noexp      = w_noexp,
    noexp_ext  = w_noexp_ext,
    notradeup  = w_notradeup,
    random     = w_random
  )
}

# ==============================================================================
# 12. MONTE CARLO: PARALLEL WITH PROGRESS BAR
# ==============================================================================

#' Create a cluster with all required exports.
#' Returned cluster should be stopped by the caller when done.
make_cluster <- function(ncores = NCORES) {
  cl <- makeCluster(ncores)
  clusterExport(cl, c(
    "init_journals", "init_articles",
    "make_signal_matrix", "build_queues", "make_crn_draws", "compute_actual_rho",
    "resolve_overbooking", "compute_review_count",
    "simulate_current", "simulate_da",
    "simulate_extended_deadlines", "simulate_no_expedite",
    "simulate_no_expedite_extended",
    "simulate_no_trading_up", "simulate_random",
    "compute_welfare", "run_one_replication",
    "STATUS_SUBMITTED", "STATUS_HOLDING", "STATUS_MATCHED"
  ))
  clusterEvalQ(cl, library(data.table))
  clusterSetRNGStream(cl, iseed = 2026)
  cl
}

run_monte_carlo <- function(p, ncores = NCORES, cl = NULL) {
  own_cluster <- is.null(cl)
  if (own_cluster) {
    message("Running ", p$n_sims, " simulations across ", ncores, " cores...")
    cl <- make_cluster(ncores)
  }

  results <- pblapply(1:p$n_sims, function(i) {
    run_one_replication(i, p)
  }, cl = cl)

  if (own_cluster) stopCluster(cl)
  message("Done.")
  results
}

# ==============================================================================
# 13. RESULTS AGGREGATION
# ==============================================================================

aggregate_results <- function(results) {
  mechanisms <- c("current", "da", "da_author", "da_constr", "extended",
                  "noexp", "noexp_ext", "notradeup", "random")
  mech_labels <- c("Current System", "Deferred Acceptance",
                    "DA (Author-Proposing)", "DA (Capacity-Constrained)",
                    "Extended Deadlines",
                    "No Expedite", "No Expedite + Extended",
                    "No Trading Up", "Random")

  # --- Aggregate match quality ---
  agg <- data.table(
    mechanism = mech_labels,
    mean_quality = sapply(mechanisms, function(m) {
      mean(sapply(results, function(r) r[[m]]$match_quality))
    }),
    sd_quality = sapply(mechanisms, function(m) {
      sd(sapply(results, function(r) r[[m]]$match_quality))
    }),
    mean_rank_cor = sapply(mechanisms, function(m) {
      mean(sapply(results, function(r) r[[m]]$rank_cor))
    }),
    mean_matched = sapply(mechanisms, function(m) {
      mean(sapply(results, function(r) r[[m]]$n_matched))
    }),
    mean_quality_allN = sapply(mechanisms, function(m) {
      mean(sapply(results, function(r) r[[m]]$match_quality_allN))
    })
  )

  # Efficiency gap relative to journal-proposing DA (all-N measure).
  # Positive = mechanism outperforms DA; negative = mechanism underperforms DA.
  da_q_allN <- agg[mechanism == "Deferred Acceptance"]$mean_quality_allN
  agg[, eff_gap_pct := round(100 * (mean_quality_allN / da_q_allN - 1), 1)]

  # --- Tier-level analysis ---
  tier_list <- list()
  for (m_idx in seq_along(mechanisms)) {
    m <- mechanisms[m_idx]
    for (tt in 1:4) {
      tier_list[[length(tier_list) + 1]] <- data.table(
        mechanism = mech_labels[m_idx],
        tier = tt,
        mean_tier_quality = mean(sapply(results, function(r) r[[m]]$tier_quality[tt])),
        mean_article_q = mean(sapply(results, function(r) r[[m]]$tier_avg_article_q[tt])),
        mean_var_q = mean(sapply(results, function(r) r[[m]]$tier_var_article_q[tt])),
        mean_n = mean(sapply(results, function(r) r[[m]]$tier_n[tt]))
      )
    }
  }
  tier_dt <- rbindlist(tier_list)

  da_tier <- tier_dt[mechanism == "Deferred Acceptance"]
  tier_dt <- merge(tier_dt, da_tier[, .(tier, da_tier_q = mean_article_q)], by = "tier")
  tier_dt[, q_gap_vs_da := mean_article_q - da_tier_q]

  # --- Author quartile analysis ---
  author_list <- list()
  for (m_idx in seq_along(mechanisms)) {
    m <- mechanisms[m_idx]
    for (qq in 1:4) {
      author_list[[length(author_list) + 1]] <- data.table(
        mechanism = mech_labels[m_idx],
        prestige_quartile = qq,
        mean_article_q = mean(sapply(results, function(r) r[[m]]$author_avg_q_by_pquartile[qq])),
        mean_j_rank = mean(sapply(results, function(r) r[[m]]$author_avg_jrank_by_pquartile[qq]))
      )
    }
  }
  author_dt <- rbindlist(author_list)

  # --- Per-tier displacement rates (fraction of tier articles that DA places
  #     in a different tier) ---
  disp_list <- list()
  for (m_idx in seq_along(mechanisms)) {
    m <- mechanisms[m_idx]
    for (tt in 1:4) {
      disp_list[[length(disp_list) + 1]] <- data.table(
        mechanism = mech_labels[m_idx],
        tier = tt,
        mean_displacement = mean(sapply(results,
          function(r) r[[m]]$displacement_vs_da[tt]))
      )
    }
  }
  displacement_dt <- rbindlist(disp_list)

  # --- Efficiency loss decomposition (spec Section 6.4) ---
  # Uses all-N quality for fair comparison (penalizes unmatched articles).
  # Total loss = DA - Current. Decompose into three additive components
  # via a telescoping path where each step changes exactly ONE feature:
  #
  #   Step 1 (expedite removal):   Current → No Expedite
  #     removes queue-jumping, same deadlines
  #   Step 2 (deadline extension):  No Expedite → No Expedite + Extended
  #     extends deadlines to 14 days, still no queue-jumping
  #   Step 3 (centralization):      No Expedite + Extended → DA
  #     replaces decentralized sequential process with DA
  #
  # The DA benchmark is journal-proposing (journal-optimal stable matching).
  q_current    <- agg[mechanism == "Current System"]$mean_quality_allN
  q_noexp      <- agg[mechanism == "No Expedite"]$mean_quality_allN
  q_noexp_ext  <- agg[mechanism == "No Expedite + Extended"]$mean_quality_allN
  q_da         <- agg[mechanism == "Deferred Acceptance"]$mean_quality_allN

  total_gap <- q_da - q_current  # negative when current system beats DA
  # Use absolute value for share computation so percentages are interpretable
  # regardless of sign. Each component's share shows what fraction of the
  # total gap it accounts for, preserving the sign of the component.
  abs_gap <- max(abs(total_gap), 1e-12)

  decomposition <- data.table(
    component = c("Expedite removal", "Deadline extension",
                  "Centralization", "Total loss vs DA (J-proposing)"),
    quality_gain = c(q_noexp - q_current,
                     q_noexp_ext - q_noexp,
                     q_da - q_noexp_ext,
                     total_gap),
    pct_of_total = round(100 * c(
      (q_noexp - q_current) / abs_gap,
      (q_noexp_ext - q_noexp) / abs_gap,
      (q_da - q_noexp_ext) / abs_gap,
      total_gap / abs_gap), 1)
  )

  # --- DA core range: difference between journal-optimal and author-optimal ---
  q_da_author <- agg[mechanism == "DA (Author-Proposing)"]$mean_quality_allN
  q_da_constr <- agg[mechanism == "DA (Capacity-Constrained)"]$mean_quality_allN
  da_core_range <- data.table(
    metric = c("Journal-Proposing DA", "Author-Proposing DA",
               "DA (Capacity-Constrained)", "Core range (J vs A)"),
    mean_quality_allN = c(q_da, q_da_author, q_da_constr,
                          abs(q_da - q_da_author))
  )

  # --- Conditional prestige-bias analysis ---
  # Aggregate the per-replication cond_prestige_bias across Monte Carlo runs.
  # For each mechanism, for each quality quintile × prestige quartile cell,
  # compute the mean journal rank across replications (ignoring NAs).
  focus_mechs <- c("current", "da")
  focus_labels <- c("Current System", "Deferred Acceptance")
  cpb_list <- list()
  for (m_idx in seq_along(focus_mechs)) {
    m <- focus_mechs[m_idx]
    # Determine number of quality bins from first replication
    n_qbins <- length(results[[1]][[m]]$cond_prestige_bias)
    for (qq in seq_len(n_qbins)) {
      for (pp in 1:4) {
        vals <- sapply(results, function(r) {
          r[[m]]$cond_prestige_bias[[qq]]$mean_jrank[pp]
        })
        ns <- sapply(results, function(r) {
          r[[m]]$cond_prestige_bias[[qq]]$n[pp]
        })
        cpb_list[[length(cpb_list) + 1]] <- data.table(
          mechanism = focus_labels[m_idx],
          quality_quintile = qq,
          prestige_quartile = pp,
          mean_jrank = mean(vals, na.rm = TRUE),
          mean_n = mean(ns, na.rm = TRUE)
        )
      }
    }
  }
  cpb_dt <- rbindlist(cpb_list)

  list(
    aggregate = agg,
    by_tier = tier_dt,
    by_author_prestige = author_dt,
    displacement = displacement_dt,
    decomposition = decomposition,
    da_core_range = da_core_range,
    cond_prestige_bias = cpb_dt
  )
}

# ==============================================================================
# 14. VISUALIZATION
# ==============================================================================

plot_results <- function(agg_results, output_dir = FIGURE_DIR) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  agg <- agg_results$aggregate
  tier_dt <- agg_results$by_tier
  author_dt <- agg_results$by_author_prestige

  tier_labels <- c("1: Elite (1-20)", "2: Top (21-50)",
                    "3: Upper-Mid (51-100)", "4: Lower (101-200)")

  five_mechs <- c("Current System", "Deferred Acceptance",
                  "No Expedite", "No Trading Up", "Random")
  five_colors <- c("grey25", "grey55", CB_PALETTE[2], CB_PALETTE[6], CB_PALETTE[7])

  # --- Fig 3: Aggregate Match Quality (dot plot) ---
  agg$mechanism <- factor(agg$mechanism, levels = agg$mechanism)

  p1 <- ggplot(agg, aes(x = mechanism, y = mean_quality_allN)) +
    geom_point(size = 3, colour = "black") +
    scale_y_continuous(limits = c(0, NA), breaks = seq(0, 0.20, 0.05)) +
    labs(x = NULL, y = "Mean Match Quality (all-N)") +
    FIG_THEME +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1, size = 9))

  save_plos_fig(file.path(output_dir, paste0("Fig3", FIG_EXT)), p1)

  # --- Fig 4: Article Quality by Journal Tier (dot plot) ---
  tier_dt$tier_label <- tier_labels[tier_dt$tier]
  tier_dt$tier_label <- factor(tier_dt$tier_label, levels = tier_labels)

  p2 <- ggplot(tier_dt[mechanism %in% c("Current System", "Deferred Acceptance")],
               aes(x = tier_label, y = mean_article_q,
                   colour = mechanism, shape = mechanism)) +
    geom_point(position = position_dodge(width = 0.4), size = 3.5) +
    scale_colour_manual(values = TWO_GROUP_COLORS) +
    scale_shape_manual(values = TWO_GROUP_SHAPES) +
    scale_y_continuous(limits = c(0, NA), breaks = seq(0, 0.6, 0.1)) +
    labs(x = "Journal Tier", y = "Mean Article Quality",
         colour = "Mechanism", shape = "Mechanism") +
    FIG_THEME +
    theme(axis.text.x = element_text(angle = 15, hjust = 1),
          legend.position = "bottom")

  save_plos_fig(file.path(output_dir, paste0("Fig4", FIG_EXT)), p2)

  # --- Fig 7: Quality Gap vs DA by Tier (dot plot) ---
  gap_mechs <- tier_dt[!mechanism %in% c("Deferred Acceptance",
                                          "DA (Author-Proposing)")]
  n_gap <- length(unique(gap_mechs$mechanism))
  gap_colors <- CB_PALETTE[1:n_gap]

  p3 <- ggplot(gap_mechs,
               aes(x = tier_label, y = q_gap_vs_da,
                   colour = mechanism, shape = mechanism)) +
    geom_point(position = position_dodge(width = 0.6), size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    scale_colour_manual(values = setNames(gap_colors, unique(gap_mechs$mechanism))) +
    scale_shape_manual(values = setNames(CB_SHAPES[1:n_gap], unique(gap_mechs$mechanism))) +
    labs(x = "Journal Tier", y = "Quality Gap (vs. DA)",
         colour = NULL, shape = NULL) +
    guides(colour = guide_legend(nrow = 4, byrow = TRUE),
           shape  = guide_legend(nrow = 4, byrow = TRUE)) +
    FIG_THEME +
    theme(axis.text.x = element_text(angle = 15, hjust = 1),
          legend.position = "bottom",
          legend.text = element_text(size = 8))

  save_plos_fig(file.path(output_dir, paste0("Fig7", FIG_EXT)), p3,
                height = FIG_H + 0.3)

  # --- Fig 9: Conditional Prestige-Bias by Quality Quintile (grouped bar) ---
  # (Was fig8 before reordering to match manuscript citation order)
  cpb <- agg_results$cond_prestige_bias
  if (!is.null(cpb) && nrow(cpb) > 0) {
    cpb[, prestige_label := factor(prestige_quartile,
      levels = 1:4,
      labels = c("Q1 (Low)", "Q2", "Q3", "Q4 (High)"))]
    cpb[, quality_label := factor(quality_quintile,
      levels = 1:5,
      labels = c("Quintile 1\n(Lowest)", "Quintile 2", "Quintile 3",
                 "Quintile 4", "Quintile 5\n(Highest)"))]

    p_cpb <- ggplot(cpb, aes(x = quality_label, y = mean_jrank,
                              fill = prestige_label)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      facet_wrap(~ mechanism, ncol = 2) +
      scale_y_reverse(name = "Mean placement rank (lower = more prestigious)") +
      scale_fill_manual(
        name = "Author\nprestige",
        values = c("Q1 (Low)" = "#333333", "Q2" = "#777777",
                   "Q3" = "#AAAAAA", "Q4 (High)" = "#DDDDDD")) +
      labs(x = "Article quality quintile") +
      FIG_THEME +
      theme(legend.position = "bottom",
            axis.text.x = element_text(size = 9),
            strip.text = element_text(size = 12, face = "bold"))

    save_plos_fig(file.path(output_dir, paste0("Fig9", FIG_EXT)), p_cpb,
                  width = 7.5, height = FIG_H + 0.5)
  } else {
    p_cpb <- NULL
  }

  # --- Fig 8: Author Prestige Quartile -> Journal Rank (dot plot) ---
  # (Was fig9 before reordering to match manuscript citation order)
  pq_labels <- c("Q1 (low prestige)", "Q2", "Q3", "Q4 (high prestige)")
  author_dt$pq_label <- pq_labels[author_dt$prestige_quartile]
  author_dt$pq_label <- factor(author_dt$pq_label, levels = pq_labels)

  p4 <- ggplot(author_dt[mechanism %in% c("Current System", "Deferred Acceptance")],
               aes(x = pq_label, y = mean_j_rank,
                   colour = mechanism, shape = mechanism)) +
    geom_point(position = position_dodge(width = 0.4), size = 3.5) +
    scale_y_reverse(limits = c(200, 1), breaks = c(1, 50, 100, 150, 200)) +
    scale_colour_manual(values = TWO_GROUP_COLORS) +
    scale_shape_manual(values = TWO_GROUP_SHAPES) +
    labs(x = "Author Prestige Quartile", y = "Mean Journal Rank (1 = most prestigious)",
         colour = "Mechanism", shape = "Mechanism") +
    FIG_THEME +
    theme(legend.position = "bottom")

  save_plos_fig(file.path(output_dir, paste0("Fig8", FIG_EXT)), p4)

  # --- Fig 5: Within-Tier Quality Variance (dot plot) ---
  p5 <- ggplot(tier_dt[mechanism %in% five_mechs],
               aes(x = tier_label, y = mean_var_q,
                   colour = mechanism, shape = mechanism)) +
    geom_point(position = position_dodge(width = 0.6), size = 3) +
    scale_colour_manual(values = setNames(five_colors, five_mechs)) +
    scale_shape_manual(values = setNames(CB_SHAPES[1:5], five_mechs)) +
    scale_y_continuous(limits = c(0, NA), breaks = seq(0, 0.030, 0.005)) +
    labs(x = "Journal Tier", y = "Within-Tier Variance\nin Article Quality",
         colour = NULL, shape = NULL) +
    guides(colour = guide_legend(nrow = 3, byrow = TRUE),
           shape  = guide_legend(nrow = 3, byrow = TRUE)) +
    FIG_THEME +
    theme(axis.text.x = element_text(angle = 15, hjust = 1),
          legend.position = "bottom",
          legend.text = element_text(size = 8))

  save_plos_fig(file.path(output_dir, paste0("Fig5", FIG_EXT)), p5,
                height = FIG_H + 0.3)

  # --- Fig 6: Displacement Rate by Tier (dot plot) ---
  disp_dt <- agg_results$displacement
  if (!is.null(disp_dt) && nrow(disp_dt) > 0) {
    disp_dt$tier_label <- tier_labels[disp_dt$tier]
    disp_dt$tier_label <- factor(disp_dt$tier_label, levels = tier_labels)

    disp_mechs <- disp_dt[!mechanism %in% c("Deferred Acceptance",
                                              "DA (Author-Proposing)")]
    n_disp <- length(unique(disp_mechs$mechanism))
    disp_colors <- CB_PALETTE[1:n_disp]

    p6 <- ggplot(disp_mechs,
                 aes(x = tier_label, y = mean_displacement,
                     colour = mechanism, shape = mechanism)) +
      geom_point(position = position_dodge(width = 0.6), size = 3) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                         labels = scales::percent_format(accuracy = 1)) +
      scale_colour_manual(values = setNames(disp_colors, unique(disp_mechs$mechanism))) +
      scale_shape_manual(values = setNames(CB_SHAPES[1:n_disp], unique(disp_mechs$mechanism))) +
      labs(x = "Journal Tier", y = "Displacement Rate",
           colour = NULL, shape = NULL) +
      guides(colour = guide_legend(nrow = 4, byrow = TRUE),
             shape  = guide_legend(nrow = 4, byrow = TRUE)) +
      FIG_THEME +
      theme(axis.text.x = element_text(angle = 15, hjust = 1),
            legend.position = "bottom",
            legend.text = element_text(size = 8))

    save_plos_fig(file.path(output_dir, paste0("Fig6", FIG_EXT)), p6,
                  height = FIG_H + 0.3)
  } else {
    p6 <- NULL
  }

  message("Figures saved to: ", output_dir)
  invisible(list(p1 = p1, p2 = p2, p3 = p3, p_cpb = p_cpb, p4 = p4, p5 = p5, p6 = p6))
}

# ==============================================================================
# 15. SENSITIVITY ANALYSIS (over rho)
# ==============================================================================

run_sensitivity_rho <- function(rho_values = seq(0.1, 0.8, by = 0.1),
                                 n_sims_each = 200,
                                 ncores = NCORES) {
  message("Sensitivity analysis over rho (mixing weight).")
  message("Mixing weights: ", paste(rho_values, collapse = ", "))
  message("Actual correlations: ",
          paste(round(compute_actual_rho(rho_values,
            params$q_alpha, params$q_beta,
            params$eps_alpha, params$eps_beta), 3), collapse = ", "))

  # Reuse a single cluster across all rho values
  cl <- make_cluster(ncores)
  on.exit(stopCluster(cl))

  all_results <- list()

  for (rv in rho_values) {
    p_mod <- params
    p_mod$rho <- rv
    p_mod$n_sims <- n_sims_each

    message("\n--- rho = ", rv,
            " (actual cor = ", round(compute_actual_rho(rv,
              p_mod$q_alpha, p_mod$q_beta,
              p_mod$eps_alpha, p_mod$eps_beta), 3), ") ---")

    mc <- run_monte_carlo(p_mod, ncores, cl = cl)
    agg <- aggregate_results(mc)
    agg$rho <- rv
    agg$actual_cor <- compute_actual_rho(rv, p_mod$q_alpha, p_mod$q_beta,
                                            p_mod$eps_alpha, p_mod$eps_beta)
    all_results[[as.character(rv)]] <- agg
  }

  # Compile key metrics across rho values
  sensitivity_dt <- rbindlist(lapply(names(all_results), function(rv) {
    a <- all_results[[rv]]$aggregate
    a$rho <- as.numeric(rv)
    a$actual_cor <- all_results[[rv]]$actual_cor
    a
  }))

  tier_sensitivity <- rbindlist(lapply(names(all_results), function(rv) {
    a <- all_results[[rv]]$by_tier
    a$rho <- as.numeric(rv)
    a$actual_cor <- all_results[[rv]]$actual_cor
    a
  }))

  out <- list(
    sensitivity_dt = sensitivity_dt,
    tier_sensitivity = tier_sensitivity,
    full_results = all_results
  )
  if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
  saveRDS(out, file.path(OUTPUT_DIR, "sensitivity_rho.rds"))
  message("Saved to ", file.path(OUTPUT_DIR, "sensitivity_rho.rds"))
  out
}

# ==============================================================================
# 15b. SENSITIVITY ANALYSIS (over aspiration scale gamma)
# ==============================================================================

run_sensitivity_gamma <- function(gamma_values = seq(0.2, 0.6, by = 0.1),
                                   n_sims_each = 200,
                                   ncores = NCORES) {
  message("Sensitivity analysis over aspiration scale (gamma).")
  message("Gamma values: ", paste(gamma_values, collapse = ", "))

  cl <- make_cluster(ncores)
  on.exit(stopCluster(cl))

  all_results <- list()

  for (gv in gamma_values) {
    p_mod <- params
    p_mod$aspiration_scale <- gv
    p_mod$n_sims <- n_sims_each

    message("\n--- aspiration_scale = ", gv, " ---")

    mc <- run_monte_carlo(p_mod, ncores, cl = cl)
    agg <- aggregate_results(mc)
    agg$aspiration_scale <- gv
    all_results[[as.character(gv)]] <- agg
  }

  sensitivity_dt <- rbindlist(lapply(names(all_results), function(gv) {
    a <- all_results[[gv]]$aggregate
    a$aspiration_scale <- as.numeric(gv)
    a
  }))

  tier_sensitivity <- rbindlist(lapply(names(all_results), function(gv) {
    a <- all_results[[gv]]$by_tier
    a$aspiration_scale <- as.numeric(gv)
    a
  }))

  out <- list(
    sensitivity_dt = sensitivity_dt,
    tier_sensitivity = tier_sensitivity,
    full_results = all_results
  )
  if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
  saveRDS(out, file.path(OUTPUT_DIR, "sensitivity_gamma.rds"))
  message("Saved to ", file.path(OUTPUT_DIR, "sensitivity_gamma.rds"))
  out
}

# ==============================================================================
# 16. ROBUSTNESS: DEADLINE BEHAVIOR
# ==============================================================================
# Runs the full Monte Carlo under both "accept" and "lapse" deadline behaviors.

run_robustness_deadline <- function(n_sims_each = 200, ncores = NCORES) {
  message("Robustness check: deadline behavior (accept vs. lapse)")

  cl <- make_cluster(ncores)
  on.exit(stopCluster(cl))

  results_out <- list()
  for (behavior in c("accept", "lapse")) {
    message("\n--- deadline_behavior = ", behavior, " ---")
    p_mod <- params
    p_mod$deadline_behavior <- behavior
    p_mod$n_sims <- n_sims_each

    mc <- run_monte_carlo(p_mod, ncores, cl = cl)
    agg <- aggregate_results(mc)
    agg$deadline_behavior <- behavior
    results_out[[behavior]] <- agg
  }

  # Compile
  comparison <- rbindlist(lapply(names(results_out), function(b) {
    a <- results_out[[b]]$aggregate
    a$deadline_behavior <- b
    a
  }))

  out <- list(comparison = comparison, full = results_out)
  if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
  saveRDS(out, file.path(OUTPUT_DIR, "robustness_deadline.rds"))
  message("Saved to ", file.path(OUTPUT_DIR, "robustness_deadline.rds"))
  out
}

# ==============================================================================
# 17. ROBUSTNESS: MARKET TIGHTNESS
# ==============================================================================
# Varies the slot_scale parameter to test sensitivity to the ratio of total
# available slots to articles. The baseline (~1545 slots for 800 articles,
# slot_scale=1.0, ratio ~1.9:1) is a relatively loose market. Tighter markets
# (lower slot_scale) create more competition and may produce different
# mechanism rankings.

run_robustness_market_tightness <- function(
    scale_values = c(0.25, 0.5, 0.75, 1.0),
    n_sims_each = 500,
    ncores = NCORES) {

  message("Robustness check: market tightness (slot_scale)")
  message("Scale values: ", paste(scale_values, collapse = ", "))

  cl <- make_cluster(ncores)
  on.exit(stopCluster(cl))

  all_results <- list()

  for (sv in scale_values) {
    p_mod <- params
    p_mod$slot_scale <- sv
    p_mod$n_sims <- n_sims_each

    # Compute approximate total slots for reporting
    # Tier 1: 20 journals, avg slots 4.5; Tier 2: 30, avg 6;
    # Tier 3: 50, avg 7.5; Tier 4: 100, avg 9. Total baseline ~1545.
    approx_slots <- round(sv * (20*4.5 + 30*6 + 50*7.5 + 100*9))
    message("\n--- slot_scale = ", sv,
            " (approx ", approx_slots, " total slots for N=", p_mod$N, ") ---")

    mc <- run_monte_carlo(p_mod, ncores, cl = cl)
    agg <- aggregate_results(mc)
    agg$slot_scale <- sv
    agg$approx_slots <- approx_slots
    all_results[[as.character(sv)]] <- agg
  }

  # Compile key metrics across scale values
  tightness_dt <- rbindlist(lapply(names(all_results), function(sv) {
    a <- all_results[[sv]]$aggregate
    a$slot_scale <- as.numeric(sv)
    a$approx_slots <- all_results[[sv]]$approx_slots
    a
  }))

  tier_tightness <- rbindlist(lapply(names(all_results), function(sv) {
    a <- all_results[[sv]]$by_tier
    a$slot_scale <- as.numeric(sv)
    a
  }))

  out <- list(
    tightness_dt = tightness_dt,
    tier_tightness = tier_tightness,
    full_results = all_results
  )
  if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
  saveRDS(out, file.path(OUTPUT_DIR, "sensitivity_market_tightness.rds"))
  message("Saved to ", file.path(OUTPUT_DIR, "sensitivity_market_tightness.rds"))
  out
}

# ==============================================================================
# 18. MAIN EXECUTION
# ==============================================================================

main <- function() {
  message(strrep("=", 60))
  message("Law Review Simulation: Starting")
  message(strrep("=", 60))

  # Report actual rho correlation
  message("rho (mixing weight) = ", params$rho,
          "  =>  actual cor(prestige, quality) = ",
          round(compute_actual_rho(params$rho, params$q_alpha, params$q_beta,
                                    params$eps_alpha, params$eps_beta), 3))

  # --- Main Monte Carlo ---
  t0 <- Sys.time()
  results <- run_monte_carlo(params)
  t1 <- Sys.time()
  message("Monte Carlo completed in ",
          round(difftime(t1, t0, units = "mins"), 1), " minutes.")

  # --- Aggregate ---
  agg <- aggregate_results(results)

  # --- Print summary ---
  cat("\n")
  message(strrep("=", 60))
  message("AGGREGATE RESULTS")
  message(strrep("=", 60))
  print(agg$aggregate)

  cat("\n")
  message("DA CORE RANGE (journal-optimal vs. author-optimal)")
  message(strrep("-", 60))
  print(agg$da_core_range)

  cat("\n")
  message("EFFICIENCY LOSS DECOMPOSITION")
  message(strrep("-", 60))
  print(agg$decomposition)

  cat("\n")
  message("TIER-LEVEL RESULTS (Current vs. DA)")
  message(strrep("-", 60))
  print(agg$by_tier[mechanism %in% c("Current System", "Deferred Acceptance")])

  cat("\n")
  message("DISPLACEMENT RATES BY TIER")
  message(strrep("-", 60))
  print(agg$displacement[mechanism %in% c("Current System", "Extended Deadlines",
                                           "No Expedite", "No Trading Up")])

  cat("\n")
  message("AUTHOR PRESTIGE QUARTILE RESULTS (Current vs. DA)")
  message(strrep("-", 60))
  print(agg$by_author_prestige[mechanism %in% c("Current System",
                                                  "Deferred Acceptance")])

  cat("\n")
  message("CONDITIONAL PRESTIGE BIAS (quality quintile x prestige quartile)")
  message(strrep("-", 60))
  message("Mean journal rank by quality quintile and prestige quartile")
  message("(Lower rank = more prestigious placement)")
  print(agg$cond_prestige_bias)

  cat("\n")
  message("DA CAPACITY-CONSTRAINED vs. UNCONSTRAINED")
  message(strrep("-", 60))
  print(agg$aggregate[mechanism %in% c("Deferred Acceptance",
                                        "DA (Capacity-Constrained)")])

  # --- Plots (output to manuscript/figures) ---
  plots <- plot_results(agg, FIGURE_DIR)

  # --- Save raw results ---
  if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
  saveRDS(results, file.path(OUTPUT_DIR, "mc_results_raw.rds"))
  saveRDS(agg, file.path(OUTPUT_DIR, "mc_results_aggregated.rds"))

  message("\nBaseline results saved.")
  invisible(list(results = results, aggregated = agg, plots = plots))
}

# ==============================================================================
# HEATMAP GRID: slot_scale × aspiration_scale
# Sweeps both parameters jointly to produce a 2D surface of DA advantage.
# ==============================================================================

run_heatmap_grid <- function(
    slot_values       = seq(0.20, 1.10, by = 0.05),
    gamma_values      = seq(0.15, 0.65, by = 0.025),
    n_sims_each       = 200,
    ncores            = NCORES) {

  message("Heatmap grid: slot_scale x aspiration_scale")
  message("  slot_scale:  ", paste(slot_values, collapse = ", "))
  message("  aspiration:  ", paste(gamma_values, collapse = ", "))
  message("  Grid cells:  ", length(slot_values) * length(gamma_values))
  message("  Reps/cell:   ", n_sims_each)

  cl <- make_cluster(ncores)
  on.exit(stopCluster(cl))

  # Pre-allocate results list
  results_list <- list()
  total_cells <- length(slot_values) * length(gamma_values)
  cell_num <- 0

  for (sv in slot_values) {
    for (gv in gamma_values) {
      cell_num <- cell_num + 1
      message(sprintf("[%d/%d] slot_scale=%.2f, aspiration_scale=%.3f",
                      cell_num, total_cells, sv, gv))

      p_mod <- params
      p_mod$slot_scale <- sv
      p_mod$aspiration_scale <- gv
      p_mod$n_sims <- n_sims_each

      mc <- run_monte_carlo(p_mod, ncores, cl = cl)
      agg <- aggregate_results(mc)

      # Extract key metrics for Current System and DA
      current <- agg$aggregate[mechanism == "Current System"]
      da      <- agg$aggregate[mechanism == "Deferred Acceptance"]

      results_list[[cell_num]] <- data.table(
        slot_scale       = sv,
        aspiration_scale = gv,
        # Current system metrics
        current_allN     = current$mean_quality_allN,
        current_quality  = current$mean_quality,
        current_rankcor  = current$mean_rank_cor,
        current_matched  = current$mean_matched,
        # DA metrics
        da_allN          = da$mean_quality_allN,
        da_quality       = da$mean_quality,
        da_rankcor       = da$mean_rank_cor,
        da_matched       = da$mean_matched,
        # Derived: DA advantage (positive = DA better)
        da_advantage_allN    = da$mean_quality_allN - current$mean_quality_allN,
        da_advantage_quality = da$mean_quality - current$mean_quality,
        da_advantage_rankcor = da$mean_rank_cor - current$mean_rank_cor,
        match_rate_gap       = current$mean_matched - da$mean_matched
      )
    }
  }

  grid_dt <- rbindlist(results_list)

  message("\nHeatmap grid complete. ", nrow(grid_dt), " cells.")
  message("DA advantage on all-N ranges from ",
          round(min(grid_dt$da_advantage_allN), 5), " to ",
          round(max(grid_dt$da_advantage_allN), 5))

  # --- Save data ---
  if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
  saveRDS(grid_dt, file.path(OUTPUT_DIR, "heatmap_grid.rds"))
  message("Grid data saved to ", file.path(OUTPUT_DIR, "heatmap_grid.rds"))

  # --- Generate heatmap figures ---
  plots <- plot_heatmap_grid(grid_dt)

  invisible(list(grid = grid_dt, plots = plots))
}


# ==============================================================================
# HEATMAP PLOTTING
# ==============================================================================

plot_heatmap_grid <- function(grid_dt) {

  if (!dir.exists(FIGURE_DIR)) dir.create(FIGURE_DIR, recursive = TRUE)
  plots <- list()

  # Shared theme for all heatmap panels
  heatmap_theme <- FIG_THEME +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.key.width = unit(0.55, "cm"),
      legend.key.height = unit(0.2, "cm"),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 7),
      legend.margin = margin(t = 0, b = 0),
      legend.box.margin = margin(t = -5),
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 10),
      plot.margin = margin(t = 5, r = 5, b = 2, l = 3)
    )

  # ---- Panel A: DA welfare advantage (all-N) with zero contour ----
  p1 <- ggplot(grid_dt, aes(x = slot_scale, y = aspiration_scale,
                             fill = da_advantage_allN)) +
    geom_tile() +
    geom_contour(aes(z = da_advantage_allN), breaks = 0,
                 colour = "black", linewidth = 0.9, linetype = "solid") +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0,
      name = expression(Delta * " Welfare"),
      breaks = c(-0.01, 0, 0.01),
      labels = c("-0.01", "0", "0.01"),
      guide = guide_colorbar(title.position = "bottom", title.hjust = 0.5)
    ) +
    labs(x = expression(lambda),
         y = expression(gamma)) +
    heatmap_theme

  plots$welfare <- p1

  # ---- Panel B: DA sorting advantage (viridis for contrast) ----
  p2 <- ggplot(grid_dt, aes(x = slot_scale, y = aspiration_scale,
                             fill = da_advantage_rankcor)) +
    geom_tile() +
    scale_fill_viridis_c(
      option = "D", direction = 1,
      name = expression(Delta * " Rank Cor."),
      breaks = pretty,
      labels = scales::label_number(accuracy = 0.01),
      guide = guide_colorbar(title.position = "bottom", title.hjust = 0.5)
    ) +
    labs(x = expression(lambda),
         y = expression(gamma)) +
    heatmap_theme

  plots$sorting <- p2

  # ---- Panel C: Match rate gap (viridis for contrast) ----
  p3 <- ggplot(grid_dt, aes(x = slot_scale, y = aspiration_scale,
                             fill = match_rate_gap)) +
    geom_tile() +
    scale_fill_viridis_c(
      option = "D", direction = 1,
      name = "Extra Matches",
      breaks = pretty,
      labels = scales::label_number(accuracy = 1),
      guide = guide_colorbar(title.position = "bottom", title.hjust = 0.5)
    ) +
    labs(x = expression(lambda),
         y = expression(gamma)) +
    heatmap_theme

  plots$matchrate <- p3

  # ---- Composite: 3-panel patchwork figure ----
  composite <- (p1 | p2 | p3) +
    plot_annotation(tag_levels = list(c("(A)", "(B)", "(C)"))) &
    theme(plot.tag = element_text(face = "bold", size = 10))

  save_plos_fig(file.path(FIGURE_DIR, paste0("Fig10", FIG_EXT)), composite,
                width = 7.5, height = 3.5)
  message("Saved composite heatmap: ", file.path(FIGURE_DIR, paste0("Fig10", FIG_EXT)))
  plots$composite <- composite

  message("All heatmap figures saved.")
  plots
}


# ==============================================================================
# CALIBRATION FIGURE (Fig 1)
# Simulates B draws of the article population and counts submissions per
# journal rank, producing a mean + 95% band plot.
# ==============================================================================

run_calibration_figure <- function(p = params, B = 500) {
  message("Calibration figure: simulating ", B, " population draws ...")

  # Matrix: rows = journal ranks 1..J, cols = replications
  sub_counts <- matrix(0L, nrow = p$J, ncol = B)

  for (b in seq_len(B)) {
    arts <- init_articles(p)
    # Count submissions: article i submits to every journal in [sub_lo, sub_hi]
    for (i in seq_len(nrow(arts))) {
      ranks <- arts$sub_lo[i]:arts$sub_hi[i]
      sub_counts[ranks, b] <- sub_counts[ranks, b] + 1L
    }
  }

  # Summarize across replications
  cal_dt <- data.table(
    rank     = 1:p$J,
    mean_sub = rowMeans(sub_counts),
    lo95     = apply(sub_counts, 1, quantile, probs = 0.025),
    hi95     = apply(sub_counts, 1, quantile, probs = 0.975)
  )

  # Georgetown reference: ~2000 submissions/year out of ~4000 articles nationally,
  # scaled to N articles in the simulation
  georgetown_rank <- 25
  georgetown_scaled <- round(2000 * (p$N / 4000))

  # --- Plot ---
  if (!dir.exists(FIGURE_DIR)) dir.create(FIGURE_DIR, recursive = TRUE)

  p1 <- ggplot(cal_dt, aes(x = rank)) +
    geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = "grey80", alpha = 0.8) +
    geom_line(aes(y = mean_sub), linewidth = 0.7, colour = "black") +
    geom_hline(yintercept = georgetown_scaled, linetype = "dashed",
               colour = "grey40", linewidth = 0.5) +
    scale_x_continuous(breaks = seq(0, 200, 50)) +
    scale_y_continuous(limits = c(0, NA), breaks = seq(0, 800, 200)) +
    labs(x = "Journal Rank", y = "Submissions Received") +
    FIG_THEME

  save_plos_fig(file.path(FIGURE_DIR, paste0("Fig1", FIG_EXT)), p1)
  message("Saved ", file.path(FIGURE_DIR, paste0("Fig1", FIG_EXT)))

  # Save data
  if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
  saveRDS(cal_dt, file.path(OUTPUT_DIR, "calibration.rds"))

  invisible(list(data = cal_dt, plot = p1))
}


# ==============================================================================
# CONVERGENCE DIAGNOSTICS
# Computes Monte Carlo Standard Errors (MCSE) for the baseline run and
# checks convergence of heatmap grid estimates at representative cells.
# ==============================================================================

run_convergence_check <- function(p = params, ncores = NCORES) {
  message("Convergence diagnostics: starting")

  if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
  if (!dir.exists(FIGURE_DIR)) dir.create(FIGURE_DIR, recursive = TRUE)

  # --- Part 1: MCSE for baseline results ---
  # Run baseline with the full n_sims and compute per-replication metrics
  message("  Part 1: Baseline MCSE (n_sims = ", p$n_sims, ")")

  cl <- make_cluster(ncores)
  on.exit(stopCluster(cl), add = TRUE)

  results <- pblapply(1:p$n_sims, function(i) {
    run_one_replication(i, p)
  }, cl = cl)

  # Extract per-replication scalars for key mechanisms
  mechanisms <- c("current", "da")
  mech_labels <- c("Current System", "Deferred Acceptance")
  metrics <- c("match_quality_allN", "rank_cor", "n_matched")
  metric_labels <- c("Welfare (all-N)", "Rank Correlation", "Articles Matched")

  mcse_list <- list()
  for (m_idx in seq_along(mechanisms)) {
    m <- mechanisms[m_idx]
    for (k in seq_along(metrics)) {
      vals <- sapply(results, function(r) r[[m]][[metrics[k]]])
      mcse_list[[length(mcse_list) + 1]] <- data.table(
        mechanism = mech_labels[m_idx],
        metric = metric_labels[k],
        mean = mean(vals),
        sd = sd(vals),
        mcse = sd(vals) / sqrt(length(vals)),
        n_sims = length(vals)
      )
    }
  }

  # Also compute MCSE for the DA advantage (paired difference)
  for (k in seq_along(metrics)) {
    vals_da <- sapply(results, function(r) r[["da"]][[metrics[k]]])
    vals_cs <- sapply(results, function(r) r[["current"]][[metrics[k]]])
    diff_vals <- vals_da - vals_cs
    mcse_list[[length(mcse_list) + 1]] <- data.table(
      mechanism = "DA - Current (paired)",
      metric = metric_labels[k],
      mean = mean(diff_vals),
      sd = sd(diff_vals),
      mcse = sd(diff_vals) / sqrt(length(diff_vals)),
      n_sims = length(diff_vals)
    )
  }

  mcse_dt <- rbindlist(mcse_list)
  mcse_dt[, mcse_pct := abs(mcse / mean) * 100]
  message("\n  MCSE summary:")
  print(mcse_dt[, .(mechanism, metric, mean = round(mean, 5),
                     mcse = round(mcse, 6), mcse_pct = round(mcse_pct, 2))])

  # --- Part 2: Heatmap convergence at representative cells ---
  # Pick cells near the contour (where precision matters most) and at corners.
  # For each cell, compute the DA welfare advantage at increasing n_sims
  # to show convergence.
  message("\n  Part 2: Heatmap convergence check")

  # Representative cells: corners + center + near-contour
  test_cells <- data.table(
    label = c("Tight/Selective", "Tight/Broad", "Loose/Selective",
              "Loose/Broad", "Center"),
    slot_scale = c(0.30, 0.30, 1.00, 1.00, 0.60),
    aspiration_scale = c(0.20, 0.60, 0.20, 0.60, 0.40)
  )

  # Run each cell at max reps and track cumulative estimates
  max_reps <- 500
  check_points <- c(50, 100, 150, 200, 300, 400, 500)

  conv_list <- list()
  for (row_i in seq_len(nrow(test_cells))) {
    tc <- test_cells[row_i]
    message("    Cell: ", tc$label,
            " (slot=", tc$slot_scale, ", gamma=", tc$aspiration_scale, ")")

    p_mod <- p
    p_mod$slot_scale <- tc$slot_scale
    p_mod$aspiration_scale <- tc$aspiration_scale
    p_mod$n_sims <- max_reps

    cell_results <- pblapply(1:max_reps, function(i) {
      run_one_replication(i, p_mod)
    }, cl = cl)

    # Extract paired DA advantage at each replication
    cell_diffs <- sapply(cell_results, function(r) {
      r[["da"]]$match_quality_allN - r[["current"]]$match_quality_allN
    })

    # Compute cumulative mean and MCSE at each checkpoint
    for (cp in check_points) {
      if (cp > max_reps) next
      sub <- cell_diffs[1:cp]
      conv_list[[length(conv_list) + 1]] <- data.table(
        label = tc$label,
        slot_scale = tc$slot_scale,
        aspiration_scale = tc$aspiration_scale,
        n_reps = cp,
        cum_mean = mean(sub),
        cum_mcse = sd(sub) / sqrt(cp)
      )
    }
  }

  conv_dt <- rbindlist(conv_list)

  # --- Convergence plot (Fig 2) ---
  pconv <- ggplot(conv_dt, aes(x = n_reps, y = cum_mean, colour = label)) +
    geom_line(linewidth = 0.7) +
    geom_ribbon(aes(ymin = cum_mean - 1.96 * cum_mcse,
                    ymax = cum_mean + 1.96 * cum_mcse,
                    fill = label), alpha = 0.12, colour = NA) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40",
               linewidth = 0.4) +
    scale_colour_manual(values = CB_PALETTE[1:5]) +
    scale_fill_manual(values = CB_PALETTE[1:5]) +
    labs(x = "Number of Replications",
         y = "DA Welfare Advantage\n(cumulative mean)",
         colour = NULL, fill = NULL) +
    guides(colour = guide_legend(nrow = 3, byrow = TRUE),
           fill   = guide_legend(nrow = 3, byrow = TRUE)) +
    FIG_THEME +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 8))

  save_plos_fig(file.path(FIGURE_DIR, paste0("Fig2", FIG_EXT)), pconv)
  message("  Saved ", file.path(FIGURE_DIR, paste0("Fig2", FIG_EXT)))

  # --- Save all diagnostics ---
  diag <- list(mcse = mcse_dt, convergence = conv_dt)
  saveRDS(diag, file.path(OUTPUT_DIR, "convergence_diagnostics.rds"))
  message("  Saved ", file.path(OUTPUT_DIR, "convergence_diagnostics.rds"))

  message("Convergence diagnostics: complete")
  invisible(diag)
}


# ==============================================================================
# 21. EVIDENCE TABLE GENERATION
# ==============================================================================
# Reads saved .rds output files and prints LaTeX for the expanded market
# tightness table (Table 2) and the robustness summary table. These are the
# same tables hardcoded in main.tex; this function regenerates them from the
# simulation output for reproducibility.

generate_evidence_tables <- function(output_dir = OUTPUT_DIR) {

  # ---- Expanded Market Tightness Table (Table 2) ----

  message("\n=== Expanded Market Tightness Table ===\n")

  tight <- readRDS(file.path(output_dir, "sensitivity_market_tightness.rds"))
  tdt <- tight$tightness_dt
  tier_tight <- tight$tier_tightness

  current_rows <- tdt[mechanism == "Current System",
    .(slot_scale, approx_slots,
      current_allN = round(mean_quality_allN, 4),
      current_rankcor = round(mean_rank_cor, 3),
      current_matched = round(mean_matched, 0))]

  da_rows <- tdt[mechanism == "Deferred Acceptance",
    .(slot_scale,
      da_allN = round(mean_quality_allN, 4),
      da_rankcor = round(mean_rank_cor, 3),
      da_matched = round(mean_matched, 0))]

  noexp_rows <- tdt[mechanism == "No Expedite",
    .(slot_scale, noexp_allN = mean_quality_allN)]
  noexp_ext_rows <- tdt[mechanism == "No Expedite + Extended",
    .(slot_scale, noexp_ext_allN = mean_quality_allN)]

  expanded_table <- merge(current_rows, da_rows, by = "slot_scale")
  expanded_table <- merge(expanded_table, noexp_rows, by = "slot_scale")
  expanded_table <- merge(expanded_table, noexp_ext_rows, by = "slot_scale")

  expanded_table[, `:=`(
    gap_pct = round(100 * (current_allN / da_allN - 1), 1),
    delta_expedite = noexp_allN - current_allN,
    delta_deadline = noexp_ext_allN - noexp_allN,
    delta_central = da_allN - noexp_ext_allN
  )]

  expanded_table[, total_gap := da_allN - current_allN]
  expanded_table[, `:=`(
    pct_expedite = round(100 * delta_expedite / ifelse(abs(total_gap) < 1e-12, 1e-12, abs(total_gap)), 1),
    pct_deadline = round(100 * delta_deadline / ifelse(abs(total_gap) < 1e-12, 1e-12, abs(total_gap)), 1),
    pct_central  = round(100 * delta_central  / ifelse(abs(total_gap) < 1e-12, 1e-12, abs(total_gap)), 1)
  )]

  tier1_current <- tier_tight[mechanism == "Current System" & tier == 1,
    .(slot_scale, t1_current_q = round(mean_article_q, 3))]
  tier1_da <- tier_tight[mechanism == "Deferred Acceptance" & tier == 1,
    .(slot_scale, t1_da_q = round(mean_article_q, 3))]
  expanded_table <- merge(expanded_table, tier1_current, by = "slot_scale")
  expanded_table <- merge(expanded_table, tier1_da, by = "slot_scale")

  disp_list <- list()
  for (sv_name in names(tight$full_results)) {
    fr <- tight$full_results[[sv_name]]
    d1 <- fr$displacement[mechanism == "Current System" & tier == 1]$mean_displacement
    disp_list[[sv_name]] <- data.table(
      slot_scale = as.numeric(sv_name),
      t1_displacement = round(100 * d1, 1)
    )
  }
  disp_dt <- rbindlist(disp_list)
  expanded_table <- merge(expanded_table, disp_dt, by = "slot_scale")

  table2_expanded <- expanded_table[order(slot_scale), .(
    slot_scale, approx_slots,
    current_allN, da_allN, gap_pct,
    current_rankcor, da_rankcor,
    current_matched, da_matched,
    pct_expedite, pct_deadline, pct_central,
    t1_current_q, t1_da_q, t1_displacement
  )]

  message("Expanded Table 2:")
  print(table2_expanded)

  # Print LaTeX
  cat("\n--- LaTeX for expanded Table 2 ---\n")
  cat("\\begin{table}[!ht]\n\\centering\n\\small\n")
  cat("\\caption{{\\bf Market tightness sensitivity (expanded).}}\n")
  cat("\\label{tab:tightness}\n")
  cat("\\begin{tabular}{@{}cccccccccc@{}}\n\\toprule\n")
  cat(" & & \\multicolumn{2}{c}{\\textbf{All-$N$ Welfare}} & & \\multicolumn{2}{c}{\\textbf{Rank Corr.}} & \\multicolumn{2}{c}{\\textbf{Tier 1}} \\\\\n")
  cat("\\cmidrule(lr){3-4} \\cmidrule(lr){6-7} \\cmidrule(lr){8-9}\n")
  cat("\\textbf{$\\lambda$} & \\textbf{Slots} & \\textbf{Current} & \\textbf{DA} & \\textbf{Gap} & \\textbf{Current} & \\textbf{DA} & \\textbf{Qual.~Gap} & \\textbf{Displ.} \\\\\n")
  cat("\\midrule\n")
  for (i in 1:nrow(table2_expanded)) {
    r <- table2_expanded[i]
    t1_gap <- round(r$t1_da_q - r$t1_current_q, 3)
    cat(sprintf("%.2f & %s & %.4f & %.4f & $%+.1f\\%%$ & %.3f & %.3f & %.3f & %.1f\\%% \\\\\n",
      r$slot_scale,
      formatC(r$approx_slots, format = "d", big.mark = "{,}"),
      r$current_allN, r$da_allN, r$gap_pct,
      r$current_rankcor, r$da_rankcor,
      t1_gap, r$t1_displacement))
  }
  cat("\\bottomrule\n\\end{tabular}\n")
  cat("\n\\medskip\n{\\footnotesize \\textbf{Decomposition (\\% of total gap):}}\n\n")
  cat("\\begin{tabular}{@{}crrr@{}}\n\\toprule\n")
  cat("\\textbf{$\\lambda$} & \\textbf{Expedite} & \\textbf{Deadline} & \\textbf{Central.} \\\\\n")
  cat("\\midrule\n")
  for (i in 1:nrow(table2_expanded)) {
    r <- table2_expanded[i]
    cat(sprintf("%.2f & %.1f\\%% & %.1f\\%% & %.1f\\%% \\\\\n",
      r$slot_scale, r$pct_expedite, r$pct_deadline, r$pct_central))
  }
  cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n")

  # ---- Robustness Summary Table ----

  message("\n=== Robustness Summary Table ===\n")

  # Panel A: rho sensitivity
  rho_sens <- readRDS(file.path(output_dir, "sensitivity_rho.rds"))
  rho_dt <- rho_sens$sensitivity_dt

  rho_current <- rho_dt[mechanism == "Current System",
    .(rho, current_allN = round(mean_quality_allN, 4),
      current_rankcor = round(mean_rank_cor, 3),
      current_matched = round(mean_matched, 0))]
  rho_da <- rho_dt[mechanism == "Deferred Acceptance",
    .(rho, da_allN = round(mean_quality_allN, 4),
      da_rankcor = round(mean_rank_cor, 3),
      da_matched = round(mean_matched, 0))]
  rho_table <- merge(rho_current, rho_da, by = "rho")
  rho_table[, gap_pct := round(100 * (current_allN / da_allN - 1), 1)]

  message("Rho sensitivity:")
  print(rho_table)

  # Panel B: deadline behavior
  deadline <- readRDS(file.path(output_dir, "robustness_deadline.rds"))
  dl_dt <- deadline$comparison

  dl_current <- dl_dt[mechanism == "Current System",
    .(deadline_behavior,
      current_allN = round(mean_quality_allN, 4),
      current_rankcor = round(mean_rank_cor, 3),
      current_matched = round(mean_matched, 0))]
  dl_da <- dl_dt[mechanism == "Deferred Acceptance",
    .(deadline_behavior,
      da_allN = round(mean_quality_allN, 4),
      da_rankcor = round(mean_rank_cor, 3),
      da_matched = round(mean_matched, 0))]
  dl_table <- merge(dl_current, dl_da, by = "deadline_behavior")
  dl_table[, gap_pct := round(100 * (current_allN / da_allN - 1), 1)]

  message("\nDeadline behavior robustness:")
  print(dl_table)

  # Print LaTeX
  cat("\n--- LaTeX for Robustness Table ---\n")
  cat("\\begin{table}[!ht]\n\\centering\n\\small\n")
  cat("\\caption{{\\bf Robustness checks.}}\n")
  cat("\\label{tab:robustness}\n")

  # Panel A
  cat("\\medskip\n{\\footnotesize \\textbf{Panel A: Prestige--quality mixing ($\\rho$)}}\n\n")
  cat("\\begin{tabular}{@{}ccccccc@{}}\n\\toprule\n")
  cat(" & \\multicolumn{2}{c}{\\textbf{All-$N$ Welfare}} & & \\multicolumn{2}{c}{\\textbf{Rank Corr.}} & \\\\\n")
  cat("\\cmidrule(lr){2-3} \\cmidrule(lr){5-6}\n")
  cat("\\textbf{$\\rho$} & \\textbf{Current} & \\textbf{DA} & \\textbf{Gap} & \\textbf{Current} & \\textbf{DA} & \\textbf{Matched} \\\\\n")
  cat("\\midrule\n")
  for (i in 1:nrow(rho_table)) {
    r <- rho_table[i]
    cat(sprintf("%.1f & %.4f & %.4f & $%+.1f\\%%$ & %.3f & %.3f & %d/%d \\\\\n",
      r$rho, r$current_allN, r$da_allN, r$gap_pct,
      r$current_rankcor, r$da_rankcor,
      r$current_matched, r$da_matched))
  }
  cat("\\bottomrule\n\\end{tabular}\n")

  # Panel B
  cat("\n\\bigskip\n{\\footnotesize \\textbf{Panel B: Deadline expiry behavior}}\n\n")
  cat("\\begin{tabular}{@{}lcccccc@{}}\n\\toprule\n")
  cat(" & \\multicolumn{2}{c}{\\textbf{All-$N$ Welfare}} & & \\multicolumn{2}{c}{\\textbf{Rank Corr.}} & \\\\\n")
  cat("\\cmidrule(lr){2-3} \\cmidrule(lr){5-6}\n")
  cat("\\textbf{Behavior} & \\textbf{Current} & \\textbf{DA} & \\textbf{Gap} & \\textbf{Current} & \\textbf{DA} & \\textbf{Matched} \\\\\n")
  cat("\\midrule\n")
  for (i in 1:nrow(dl_table)) {
    r <- dl_table[i]
    cat(sprintf("%s & %.4f & %.4f & $%+.1f\\%%$ & %.3f & %.3f & %d/%d \\\\\n",
      tools::toTitleCase(r$deadline_behavior),
      r$current_allN, r$da_allN, r$gap_pct,
      r$current_rankcor, r$da_rankcor,
      r$current_matched, r$da_matched))
  }
  cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n")

  message("\n=== All evidence tables generated. ===")
  message("Copy LaTeX output above into manuscript as needed.")

  invisible(list(tightness = table2_expanded, rho = rho_table, deadline = dl_table))
}


# ==============================================================================
# FULL PIPELINE EXECUTION
# Sourcing this file runs everything: baseline Monte Carlo, all sensitivity
# analyses, robustness checks, and the heatmap grid.
# All RDS outputs go to output/, all figures to manuscript/figures/.
# ==============================================================================

# Set RUN_PIPELINE <- TRUE before sourcing to run the full pipeline.
# Set RUN_PIPELINE <- "baseline" to run only the baseline Monte Carlo.
# Default: do NOT run anything (just load function definitions).
if (!exists("RUN_PIPELINE")) RUN_PIPELINE <- FALSE

if (identical(RUN_PIPELINE, TRUE) || identical(RUN_PIPELINE, "baseline")) {

  message("\n", strrep("=", 60))
  message(if (identical(RUN_PIPELINE, TRUE)) "FULL PIPELINE: Starting"
          else "BASELINE ONLY: Starting")
  message(strrep("=", 60))

  t_start <- Sys.time()

  # 0. Calibration figure (Fig1.tif, calibration.rds)
  message("\n", strrep("=", 60))
  message("CALIBRATION FIGURE")
  message(strrep("=", 60))
  calibration <- run_calibration_figure()

  # 1. Baseline Monte Carlo (produces mc_results_raw.rds, mc_results_aggregated.rds,
  #    and figures: Fig3-Fig9.tif)
  baseline <- main()

  if (identical(RUN_PIPELINE, TRUE)) {

    # 2. Sensitivity: prestige-quality correlation (rho)
    #    Produces sensitivity_rho.rds
    message("\n", strrep("=", 60))
    message("SENSITIVITY: Prestige-quality correlation (rho)")
    message(strrep("=", 60))
    sens_rho <- run_sensitivity_rho()

    # 3. Sensitivity: aspiration breadth (gamma)
    #    Produces sensitivity_gamma.rds
    message("\n", strrep("=", 60))
    message("SENSITIVITY: Aspiration breadth (gamma)")
    message(strrep("=", 60))
    sens_gamma <- run_sensitivity_gamma()

    # 4. Robustness: deadline behavior
    #    Produces robustness_deadline.rds
    message("\n", strrep("=", 60))
    message("ROBUSTNESS: Deadline behavior")
    message(strrep("=", 60))
    robust_deadline <- run_robustness_deadline()

    # 5. Robustness: market tightness (slot_scale)
    #    Produces sensitivity_market_tightness.rds
    message("\n", strrep("=", 60))
    message("ROBUSTNESS: Market tightness (slot_scale)")
    message(strrep("=", 60))
    robust_tightness <- run_robustness_market_tightness()

    # 6. Heatmap grid: slot_scale × aspiration_scale
    #    Produces heatmap_grid.rds and Fig10.tif (3-panel patchwork composite)
    message("\n", strrep("=", 60))
    message("HEATMAP GRID: slot_scale x aspiration_scale")
    message(strrep("=", 60))
    heatmap <- run_heatmap_grid()

    # 7. Convergence diagnostics (Fig2.tif, convergence_diagnostics.rds)
    #    MCSE table for baseline + convergence plot for representative heatmap cells
    message("\n", strrep("=", 60))
    message("CONVERGENCE DIAGNOSTICS")
    message(strrep("=", 60))
    convergence <- run_convergence_check()

    # 8. Evidence tables (LaTeX for expanded Table 2 and robustness table)
    #    Reads saved .rds files and prints LaTeX to console
    message("\n", strrep("=", 60))
    message("EVIDENCE TABLES (LaTeX output)")
    message(strrep("=", 60))
    evidence_tables <- generate_evidence_tables()

  }  # end full pipeline

  t_end <- Sys.time()
  message("\n", strrep("=", 60))
  message("Complete. Total time: ",
          round(difftime(t_end, t_start, units = "mins"), 1), " minutes")
  message("Output files: ", paste(list.files(OUTPUT_DIR), collapse = ", "))
  message("Figures: ", paste(list.files(FIGURE_DIR), collapse = ", "))
  message(strrep("=", 60))

}  # end RUN_PIPELINE guard
