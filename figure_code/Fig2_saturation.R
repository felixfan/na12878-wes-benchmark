# Figure 2 - Depth-saturation of variant-calling accuracy (DeepVariant, GRCh38)
# Reads source_data/Fig2_saturation.csv and writes figures/Fig2_saturation.png (300 dpi).
#   Rscript figure_code/Fig2_saturation.R
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

PLAT <- c(illumina_novaseq6000 = "#4C72B0", bgi_dnbseq_t7 = "#DD8452", genemind_surfseq5000 = "#55A868")

d <- read.csv(src("Fig2_saturation.csv"))
d <- d[d$caller == "deepvariant", ]
d$variant_type <- factor(d$variant_type, levels = c("SNV", "INDEL"))

p <- ggplot(d, aes(mean_on_target_depth, F1_mean, colour = platform)) +
  geom_line(linewidth = .7) +
  geom_point(size = 1.7) +
  geom_errorbar(aes(ymin = F1_mean - F1_std, ymax = F1_mean + F1_std), width = .03, na.rm = TRUE) +
  facet_wrap(~variant_type, scales = "free_y") +
  scale_x_log10(breaks = c(8, 10, 15, 20, 30, 50, 75, 100, 200, 400, 700)) +
  scale_colour_manual(values = PLAT) +
  labs(title = "Figure 2 - Depth-saturation of variant-calling accuracy (DeepVariant)",
       subtitle = "F1 versus measured on-target depth; error bars are +/- 1 SD over three downsampling seeds",
       x = "Measured on-target mean depth (X, log scale)", y = "F1", colour = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(out("Fig2_saturation.png"), p, width = 11, height = 4.8, dpi = 300)
cat("->", out("Fig2_saturation.png"), "\n")
