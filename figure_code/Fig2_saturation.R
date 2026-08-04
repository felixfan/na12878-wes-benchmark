# Figure 2 - Depth saturation of variant-calling accuracy (DeepVariant, GRCh38)
# Reads source_data/Fig2_saturation.csv and writes figures/fig2.png (300 dpi).
#   Rscript figure_code/Fig2_saturation.R
# Panels: (a) SNV, (b) INDEL. No title is drawn on the figure - it belongs in the legend.
# Runs from any working directory. Needs: ggplot2, patchwork
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})

.here <- local({
  a <- commandArgs(FALSE); m <- grep("^--file=", a)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", a[m[1]]))) else getwd()
})
.root <- local({
  for (p in c(file.path(.here, ".."), .here, getwd()))
    if (dir.exists(file.path(p, "source_data"))) return(normalizePath(p))
  stop("cannot locate source_data/ - run this script from inside the repository")
})
src <- function(f) file.path(.root, "source_data", f)
out <- function(f) { d <- file.path(.root, "figures")
                     dir.create(d, showWarnings = FALSE, recursive = TRUE); file.path(d, f) }

PLAT_KEY <- c("illumina_novaseq6000", "bgi_dnbseq_t7", "genemind_surfseq5000")
PLAT_LAB <- c("Illumina NovaSeq 6000", "BGI DNBSEQ-T7", "GeneMind SURFSeq 5000")
PLAT_COL <- setNames(c("#4C72B0", "#DD8452", "#55A868"), PLAT_LAB)

d <- read.csv(src("Fig2_saturation.csv"))
d <- d[d$caller == "deepvariant", ]
d$platform <- factor(PLAT_LAB[match(d$platform, PLAT_KEY)], levels = PLAT_LAB)

panel <- function(vt) {
  ggplot(d[d$variant_type == vt, ], aes(mean_on_target_depth, F1_mean, colour = platform)) +
    geom_line(linewidth = .7) +
    geom_point(size = 1.7) +
    geom_errorbar(aes(ymin = F1_mean - F1_std, ymax = F1_mean + F1_std), width = .03, na.rm = TRUE) +
    scale_x_log10(breaks = c(8, 10, 15, 20, 30, 50, 75, 100, 200, 400, 700)) +
    scale_colour_manual(values = PLAT_COL) +
    labs(x = "Measured on-target mean depth (X, log scale)", y = "F1", colour = NULL) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
}

p <- (panel("SNV") | panel("INDEL")) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(legend.position = "bottom", plot.tag = element_text(face = "bold", size = 12))

ggsave(out("fig2.png"), p, width = 11, height = 4.8, dpi = 300)
cat("->", out("fig2.png"), "\n")
