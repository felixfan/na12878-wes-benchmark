# Figure 3 - Cross-platform comparison at nominal 30x (DeepVariant, GRCh38)
# Reads source_data/Fig3_platform.csv and writes figures/fig3.png (300 dpi).
#   Rscript figure_code/Fig3_platform.R
# Panels: (a) precision, (b) recall, (c) F1. No title is drawn on the figure.
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

d <- read.csv(src("Fig3_platform.csv"))
d <- d[d$caller == "deepvariant", ]
d$platform <- factor(PLAT_LAB[match(d$platform, PLAT_KEY)], levels = PLAT_LAB)
d$variant_type <- factor(d$variant_type, levels = c("SNV", "INDEL"))

panel <- function(metric, ylab) {
  d$value <- d[[metric]]
  ggplot(d, aes(variant_type, value, fill = platform)) +
    geom_col(position = position_dodge(.8), width = .7) +
    geom_text(aes(label = sprintf("%.3f", value)), position = position_dodge(.8),
              vjust = -0.4, size = 2.5) +
    scale_fill_manual(values = PLAT_COL) +
    coord_cartesian(ylim = c(0.8, 1.02)) +
    labs(x = NULL, y = ylab, fill = NULL) +
    theme_bw(base_size = 11) +
    theme(panel.grid.major.x = element_blank())
}

p <- (panel("precision", "Precision") | panel("recall", "Recall") | panel("F1", "F1")) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(legend.position = "bottom", plot.tag = element_text(face = "bold", size = 12))

ggsave(out("fig3.png"), p, width = 11, height = 4.4, dpi = 300)
cat("->", out("fig3.png"), "\n")
