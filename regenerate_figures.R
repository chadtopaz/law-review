# ==============================================================================
# regenerate_figures.R
#
# Regenerates all 10 manuscript figures from saved simulation output (.rds files)
# WITHOUT re-running any simulations. Produces PLOS ONE-compliant TIFF files.
#
# Output: manuscript/plos/figures/Fig1.tif through Fig10.tif
#
# NOTE: Fig 8 and Fig 9 were swapped relative to original numbering to match
# manuscript citation order (author prestige = Fig 8, cond. prestige = Fig 9).
#
# Usage:
#   setwd("<project root>")   # directory containing output/ and manuscript/
#   source("code/regenerate_figures.R")
#
# Prerequisites:
#   - Saved .rds files in output/ (produced by simulation.R)
#   - R packages: data.table, ggplot2, scales, patchwork
# ==============================================================================

library(data.table)
library(ggplot2)
library(scales)
library(patchwork)

# --- Paths ---
OUTPUT_DIR <- "output"
FIGURE_DIR <- file.path("manuscript", "plos", "figures")
if (!dir.exists(FIGURE_DIR)) dir.create(FIGURE_DIR, recursive = TRUE)

# --- PLOS ONE figure specifications ---
FIG_W   <- 5.2       # column width in inches
FIG_H   <- 3.6       # default height in inches
FIG_DPI <- 300        # minimum 300, maximum 600
FIG_FORMAT <- "tiff"
FIG_EXT <- ".tif"

# --- Styling constants ---
CB_PALETTE <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
CB_SHAPES <- c(16, 17, 15, 18, 8, 4, 3)

TWO_GROUP_COLORS <- c("Current System" = "grey25", "Deferred Acceptance" = "grey60")
TWO_GROUP_SHAPES <- c("Current System" = 16, "Deferred Acceptance" = 17)

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
  message("  Saved: ", filename)
}

# ==============================================================================
# Load saved data
# ==============================================================================
message("Loading saved simulation results ...")

cal_dt   <- readRDS(file.path(OUTPUT_DIR, "calibration.rds"))
diag     <- readRDS(file.path(OUTPUT_DIR, "convergence_diagnostics.rds"))
agg_all  <- readRDS(file.path(OUTPUT_DIR, "mc_results_aggregated.rds"))
grid_dt  <- readRDS(file.path(OUTPUT_DIR, "heatmap_grid.rds"))

agg       <- agg_all$aggregate
tier_dt   <- agg_all$by_tier
author_dt <- agg_all$by_author_prestige
cpb       <- agg_all$cond_prestige_bias
disp_dt   <- agg_all$displacement
conv_dt   <- diag$convergence

message("  All data loaded.")

# ==============================================================================
# Shared labels and palettes
# ==============================================================================
tier_labels <- c("1: Elite (1-20)", "2: Top (21-50)",
                  "3: Upper-Mid (51-100)", "4: Lower (101-200)")
tier_dt$tier_label <- factor(tier_labels[tier_dt$tier], levels = tier_labels)

five_mechs  <- c("Current System", "Deferred Acceptance",
                 "No Expedite", "No Trading Up", "Random")
five_colors <- c("grey25", "grey55", CB_PALETTE[2], CB_PALETTE[6], CB_PALETTE[7])

# ==============================================================================
# Fig 1: Calibration — Submissions per journal rank
# ==============================================================================
message("Fig 1: Calibration ...")
p_cal <- ggplot(cal_dt, aes(x = rank)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = "grey80", alpha = 0.8) +
  geom_line(aes(y = mean_sub), linewidth = 0.7, colour = "black") +
  geom_hline(yintercept = 400, linetype = "dashed",
             colour = "grey40", linewidth = 0.5) +
  scale_x_continuous(breaks = seq(0, 200, 50)) +
  scale_y_continuous(limits = c(0, NA), breaks = seq(0, 800, 200)) +
  labs(x = "Journal Rank", y = "Submissions Received") +
  FIG_THEME

save_plos_fig(file.path(FIGURE_DIR, paste0("Fig1", FIG_EXT)), p_cal)

# ==============================================================================
# Fig 2: Convergence diagnostics
# ==============================================================================
message("Fig 2: Convergence ...")
p_conv <- ggplot(conv_dt, aes(x = n_reps, y = cum_mean, colour = label)) +
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

save_plos_fig(file.path(FIGURE_DIR, paste0("Fig2", FIG_EXT)), p_conv)

# ==============================================================================
# Fig 3: Aggregate match quality by mechanism (dot plot)
# ==============================================================================
message("Fig 3: Aggregate match quality ...")
agg$mechanism <- factor(agg$mechanism, levels = agg$mechanism)

p_agg <- ggplot(agg, aes(x = mechanism, y = mean_quality_allN)) +
  geom_point(size = 3, colour = "black") +
  scale_y_continuous(limits = c(0, NA), breaks = seq(0, 0.20, 0.05)) +
  labs(x = NULL, y = "Mean Match Quality (all-N)") +
  FIG_THEME +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1, size = 9))

save_plos_fig(file.path(FIGURE_DIR, paste0("Fig3", FIG_EXT)), p_agg)

# ==============================================================================
# Fig 4: Article quality by journal tier
# ==============================================================================
message("Fig 4: Quality by tier ...")
p_tier <- ggplot(tier_dt[mechanism %in% c("Current System", "Deferred Acceptance")],
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

save_plos_fig(file.path(FIGURE_DIR, paste0("Fig4", FIG_EXT)), p_tier)

# ==============================================================================
# Fig 5: Within-tier quality variance
# ==============================================================================
message("Fig 5: Within-tier variance ...")
p_var <- ggplot(tier_dt[mechanism %in% five_mechs],
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

save_plos_fig(file.path(FIGURE_DIR, paste0("Fig5", FIG_EXT)), p_var,
              height = FIG_H + 0.3)

# ==============================================================================
# Fig 6: Displacement rate by tier
# ==============================================================================
message("Fig 6: Displacement rates ...")
if (!is.null(disp_dt) && nrow(disp_dt) > 0) {
  disp_dt$tier_label <- factor(tier_labels[disp_dt$tier], levels = tier_labels)
  disp_mechs <- disp_dt[!mechanism %in% c("Deferred Acceptance",
                                            "DA (Author-Proposing)")]
  n_disp <- length(unique(disp_mechs$mechanism))
  disp_colors <- CB_PALETTE[1:n_disp]

  p_disp <- ggplot(disp_mechs,
                   aes(x = tier_label, y = mean_displacement,
                       colour = mechanism, shape = mechanism)) +
    geom_point(position = position_dodge(width = 0.6), size = 3) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                       labels = percent_format(accuracy = 1)) +
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

  save_plos_fig(file.path(FIGURE_DIR, paste0("Fig6", FIG_EXT)), p_disp,
                height = FIG_H + 0.3)
} else {
  message("  Skipped — no displacement data.")
}

# ==============================================================================
# Fig 7: Quality gap vs DA by tier
# ==============================================================================
message("Fig 7: Quality gap ...")
gap_mechs <- tier_dt[!mechanism %in% c("Deferred Acceptance",
                                        "DA (Author-Proposing)")]
n_gap <- length(unique(gap_mechs$mechanism))
gap_colors <- CB_PALETTE[1:n_gap]

p_gap <- ggplot(gap_mechs,
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

save_plos_fig(file.path(FIGURE_DIR, paste0("Fig7", FIG_EXT)), p_gap,
              height = FIG_H + 0.3)

# ==============================================================================
# Fig 8: Mean placement rank by author prestige quartile (dot plot)
# (Cited first in manuscript; was formerly fig9)
# ==============================================================================
message("Fig 8: Author prestige placement ...")
pq_labels <- c("Q1 (low prestige)", "Q2", "Q3", "Q4 (high prestige)")
author_dt$pq_label <- factor(pq_labels[author_dt$prestige_quartile],
                              levels = pq_labels)

p_auth <- ggplot(author_dt[mechanism %in% c("Current System", "Deferred Acceptance")],
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

save_plos_fig(file.path(FIGURE_DIR, paste0("Fig8", FIG_EXT)), p_auth)

# ==============================================================================
# Fig 9: Conditional prestige bias by quality quintile (grouped bar)
# (Cited second in manuscript; was formerly fig8)
# ==============================================================================
message("Fig 9: Conditional prestige bias ...")
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
          axis.text.x = element_text(size = 8),
          strip.text = element_text(size = 10, face = "bold"))

  save_plos_fig(file.path(FIGURE_DIR, paste0("Fig9", FIG_EXT)), p_cpb,
                width = 7.5, height = FIG_H + 0.5)
} else {
  message("  Skipped — no conditional prestige data.")
}

# ==============================================================================
# Fig 10: Heatmap composite (welfare gap, rank correlation gap, extra matches)
# ==============================================================================
message("Fig 10: Heatmap composite ...")

heatmap_theme <- FIG_THEME +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        legend.key.width = unit(0.55, "cm"),
        legend.key.height = unit(0.2, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.margin = margin(t = 0, b = 0),
        legend.box.margin = margin(t = -5),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 10),
        plot.margin = margin(t = 5, r = 5, b = 2, l = 3))

# Panel A: DA welfare advantage (all-N) with zero contour
p_heat1 <- ggplot(grid_dt, aes(x = slot_scale, y = aspiration_scale,
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
  labs(x = expression(lambda), y = expression(gamma)) +
  heatmap_theme

# Panel B: DA sorting advantage
p_heat2 <- ggplot(grid_dt, aes(x = slot_scale, y = aspiration_scale,
                                fill = da_advantage_rankcor)) +
  geom_tile() +
  scale_fill_viridis_c(
    option = "D", direction = 1,
    name = expression(Delta * " Rank Cor."),
    breaks = pretty,
    labels = label_number(accuracy = 0.01),
    guide = guide_colorbar(title.position = "bottom", title.hjust = 0.5)
  ) +
  labs(x = expression(lambda), y = expression(gamma)) +
  heatmap_theme

# Panel C: Match rate gap
p_heat3 <- ggplot(grid_dt, aes(x = slot_scale, y = aspiration_scale,
                                fill = match_rate_gap)) +
  geom_tile() +
  scale_fill_viridis_c(
    option = "D", direction = 1,
    name = "Extra Matches",
    breaks = pretty,
    labels = label_number(accuracy = 1),
    guide = guide_colorbar(title.position = "bottom", title.hjust = 0.5)
  ) +
  labs(x = expression(lambda), y = expression(gamma)) +
  heatmap_theme

composite <- (p_heat1 | p_heat2 | p_heat3) +
  plot_annotation(tag_levels = list(c("(A)", "(B)", "(C)"))) &
  theme(plot.tag = element_text(face = "bold", size = 10))

save_plos_fig(file.path(FIGURE_DIR, paste0("Fig10", FIG_EXT)), composite,
              width = 7.5, height = 3.5)

# ==============================================================================
message("\nAll 10 figures regenerated in: ", FIGURE_DIR)
message("Files: Fig1.tif through Fig10.tif")
