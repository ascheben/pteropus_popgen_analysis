# Compare observed and fastsimcoal2-simulated (expected) joint 2D site frequency
# spectra (SFS) for the best replicate of the best-fitting model per population pair
# (human mutation rate).
#
# Run from this directory:  cd 5_Fastsimcoal/scripts && Rscript obs_vs_exp.R
#
# Inputs:
#   ../sfs/<pair>_jointMAFpop1_0.obs                     observed SFS (easySFS)
#   ../data/human_rate/best_replicates/<rep>/<rep>_jointMAFpop1_0.txt
#                                                        expected SFS (fsc --initValues ... OUTEXP)
# Outputs (panels; titles, axis labels, panel letters and legend were finalised
# in a graphics editor for the manuscript):
#   ../figures/2D_SFS_basic_plot_small.pdf   Figure S17 (Table 1 population pairs)
#   ../figures/2D_SFS_plot_big_poor_fit.pdf  Figure S18 (all Australian P. alecto vs Indonesia)

library(ggplot2)
library(reshape2)
library(cowplot)

safe_log10 <- function(x) {
  x <- as.matrix(x)
  x[!is.finite(x)] <- 0
  x[x < 0] <- 0
  log10(x + 1)
}

# Read observed and expected SFS, scale expected to the observed number of
# polymorphic sites and return long-format data frames (log10 counts).
read_sfs_pair <- function(exp_file, obs_file) {

  obs <- as.matrix(
    read.table(obs_file, skip = 1, header = TRUE, check.names = FALSE,
               comment.char = "", stringsAsFactors = FALSE)
  )
  exp <- as.matrix(
    read.table(exp_file, header = TRUE, check.names = FALSE,
               comment.char = "", stringsAsFactors = FALSE)
  )

  # Normalize expected SFS to observed total
  obs_total <- sum(obs, na.rm = TRUE) - obs["d1_0","d0_0"]
  exp_total <- sum(exp, na.rm = TRUE)
  exp <- exp * obs_total / exp_total

  to_df <- function(m) {
    df <- melt(safe_log10(m))
    colnames(df) <- c("Pop1", "Pop2", "Value")
    # Extract numeric SFS coordinates
    df$Pop1 <- as.numeric(gsub("d1_", "", df$Pop1))
    df$Pop2 <- as.numeric(gsub("d0_", "", df$Pop2))
    # Remove unestimated ancestral class
    df[df$Pop1 == 0 & df$Pop2 == 0,] <- NA
    df
  }

  list(obs = to_df(obs), exp = to_df(exp))
}

# Tile plot of the first 15 x 15 SFS classes
sfs_tile <- function(df) {
  ggplot(
    df[df$Pop1 < 15 & df$Pop2 < 15,],
    aes(x = Pop1, y = Pop2, fill = Value)
  ) +
    geom_tile() +
    coord_equal() +
    scale_fill_viridis_b(
      name = expression(log[10](count + 1)),
      breaks = seq(0, 6, by = 1)
    ) +
    theme_bw()
}

# Observed vs expected with titles and legend
make_sfs_plots <- function(exp_file, obs_file, title_prefix = "") {
  sfs <- read_sfs_pair(exp_file, obs_file)
  style <- function(g, what) {
    g +
      theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
      labs(
        title = paste(title_prefix, what),
        x = "SNP count population 1",
        y = "SNP count population 2"
      )
  }
  plot_grid(style(sfs_tile(sfs$obs), "Observed 2D-SFS"),
            style(sfs_tile(sfs$exp), "Expected 2D-SFS"))
}

# Observed vs expected without titles, axis labels or legend (compact panels)
make_basic_sfs_plots <- function(exp_file, obs_file) {
  sfs <- read_sfs_pair(exp_file, obs_file)
  style <- function(g) {
    g +
      theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt"),
            legend.position = "none",
            plot.title = element_text(hjust = 0.5, face = "bold")) +
      labs(x = "", y = "")
  }
  plot_grid(style(sfs_tile(sfs$obs)), style(sfs_tile(sfs$exp)))
}

# Locate the observed and expected SFS for the best replicate of a pair
best_rep_dir <- "../data/human_rate/best_replicates"
sfs_files <- function(pair) {
  rep <- list.files(best_rep_dir, pattern = paste0("^", pair, "_m[1-4]_.*_rep[0-9]+$"))
  stopifnot(length(rep) == 1)
  list(
    exp = file.path(best_rep_dir, rep, paste0(rep, "_jointMAFpop1_0.txt")),
    obs = file.path("..", "sfs", paste0(pair, "_jointMAFpop1_0.obs")),
    name = pair
  )
}

dir.create("../figures", showWarnings = FALSE)

# Figure S17: population pairs reported in Table 1
# Pop codes: BFFINDOplusWAmuseum = INDO+3NWA, BFFNAUS = NWA+BURK, BFFEAST = EA+NEA, SFF = WTA+NG
table1_pairs <- c("BFFEAST_SFF",
                  "BFFNAUS_SFF",
                  "BFFINDOplusWAmuseum_BFFNAUS",
                  "BFFINDOplusWAmuseum_BFFEAST",
                  "BFFNAUS_BFFEAST",
                  "BFFINDOplusWAmuseum_SFF")
basic_plots_table1 <- lapply(table1_pairs, function(p) {
  f <- sfs_files(p)
  make_basic_sfs_plots(exp_file = f$exp, obs_file = f$obs)
})
plot_grid(plotlist = basic_plots_table1, ncol = 2)
ggsave("../figures/2D_SFS_basic_plot_small.pdf", height = 20, width = 14, dpi = 600)

# Figure S18: all Australian P. alecto (BFFAUS = NWA+EA+NEA) vs Indonesia
# A) with (BFFINDOplusWAmuseum = INDO+3NWA) and B) without (BFFINDO = INDO)
# the three North West Australian samples clustering with Indonesia
poor_fit_plots <- lapply(c("BFFAUS_BFFINDOplusWAmuseum", "BFFAUS_BFFINDO"), function(p) {
  f <- sfs_files(p)
  make_sfs_plots(exp_file = f$exp, obs_file = f$obs, title_prefix = f$name)
})
plot_grid(plotlist = poor_fit_plots, ncol = 1)
ggsave("../figures/2D_SFS_plot_big_poor_fit.pdf", height = 8, width = 8, dpi = 300)
