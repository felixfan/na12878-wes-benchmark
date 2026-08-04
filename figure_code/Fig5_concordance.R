# Figure 5 - CPU versus GPU concordance, 132 matched pairs (GRCh38)
# Reads source_data/Fig5_concordance.csv and writes figures/Fig5_concordance.png (300 dpi).
#   Rscript figure_code/Fig5_concordance.R
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

d <- read.csv(src("Fig5_concordance.csv"))
lo <- min(c(d$CPU_F1, d$GPU_F1)) - 0.004
CALLER <- c(gatk4_hc = "#4C72B0", deepvariant = "#DD8452")

p <- ggplot(d, aes(CPU_F1, GPU_F1, colour = caller_pair, shape = variant_type)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50", linewidth = .4) +
  geom_point(size = 2, alpha = .85) +
  scale_colour_manual(values = CALLER,
                      labels = c(deepvariant = "DeepVariant", gatk4_hc = "HaplotypeCaller")) +
  coord_equal(xlim = c(lo, 1.001), ylim = c(lo, 1.001)) +
  labs(title = "Figure 5 - CPU versus GPU concordance (F1)",
       subtitle = "132 matched pairs, SNV and INDEL; the dashed line is y = x",
       x = "CPU F1", y = "GPU F1 (NVIDIA Parabricks)", colour = NULL, shape = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right", panel.grid.minor = element_blank())

ggsave(out("Fig5_concordance.png"), p, width = 7, height = 5.4, dpi = 300)
cat("->", out("Fig5_concordance.png"), "\n")
