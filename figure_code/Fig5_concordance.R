# Figure 5 - CPU versus GPU concordance, 132 matched pairs (GRCh38)
# Reads source_data/Fig5_concordance.csv and writes figures/fig5.png (300 dpi).
#   Rscript figure_code/Fig5_concordance.R
# Single panel. No title is drawn on the figure - it belongs in the legend.
# Runs from any working directory. Needs: ggplot2
suppressPackageStartupMessages(library(ggplot2))

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

CALLER_LAB <- c(gatk4_hc = "GATK4 HaplotypeCaller", deepvariant = "DeepVariant")
CALLER_COL <- c("GATK4 HaplotypeCaller" = "#4C72B0", "DeepVariant" = "#DD8452")

d <- read.csv(src("Fig5_concordance.csv"))
d$caller_pair <- factor(CALLER_LAB[d$caller_pair], levels = names(CALLER_COL))
d$variant_type <- factor(d$variant_type, levels = c("SNV", "INDEL"))
lo <- min(c(d$CPU_F1, d$GPU_F1)) - 0.004

p <- ggplot(d, aes(CPU_F1, GPU_F1, colour = caller_pair, shape = variant_type)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50", linewidth = .4) +
  geom_point(size = 2, alpha = .85) +
  scale_colour_manual(values = CALLER_COL) +
  coord_equal(xlim = c(lo, 1.001), ylim = c(lo, 1.001)) +
  labs(x = "CPU F1", y = "GPU F1 (NVIDIA Parabricks)", colour = NULL, shape = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right", panel.grid.minor = element_blank())

ggsave(out("fig5.png"), p, width = 7, height = 5.2, dpi = 300)
cat("->", out("fig5.png"), "\n")
