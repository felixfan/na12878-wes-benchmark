# Figure 4 - GATK INDEL-only hard-filtering, before and after (GATK4 HaplotypeCaller, full depth, GRCh38)
# Reads source_data/Fig4_hardfilter.csv and writes figures/fig4.png (300 dpi).
#   Rscript figure_code/Fig4_hardfilter.R
# Panels: (a) F1, (b) precision, (c) recall. No title is drawn on the figure.
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

FILT_COL <- c("Raw (unfiltered)" = "#B0B0B0", "INDEL hard-filtered" = "#4C72B0")

d <- read.csv(src("Fig4_hardfilter.csv"), check.names = FALSE)
d$variant_type <- factor(d$variant_type, levels = c("SNV", "INDEL"))
d$filter <- factor(c("raw (unfiltered)" = "Raw (unfiltered)",
                     "hard-filtered" = "INDEL hard-filtered")[d$filter],
                   levels = names(FILT_COL))

panel <- function(metric, ylab) {
  d$value <- d[[metric]]
  ggplot(d, aes(variant_type, value, fill = filter)) +
    geom_col(position = position_dodge(.75), width = .65) +
    geom_text(aes(label = sprintf("%.3f", value)), position = position_dodge(.75),
              vjust = -0.4, size = 2.6) +
    scale_fill_manual(values = FILT_COL) +
    coord_cartesian(ylim = c(0.6, 1.05)) +
    labs(x = NULL, y = ylab, fill = NULL) +
    theme_bw(base_size = 11) +
    theme(panel.grid.major.x = element_blank())
}

p <- (panel("F1", "F1") | panel("precision", "Precision") | panel("recall", "Recall")) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(legend.position = "bottom", plot.tag = element_text(face = "bold", size = 12))

ggsave(out("fig4.png"), p, width = 11, height = 4.4, dpi = 300)
cat("->", out("fig4.png"), "\n")
