# Figure 3 - Cross-platform comparison at nominal 30x (DeepVariant, GRCh38)
# Reads source_data/Fig3_platform.csv and writes figures/Fig3_platform.png (300 dpi).
#   Rscript figure_code/Fig3_platform.R
# Runs from any working directory. Needs: ggplot2, tidyr
suppressPackageStartupMessages({library(ggplot2); library(tidyr)})

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

d <- read.csv(src("Fig3_platform.csv"))
d <- d[d$caller == "deepvariant", ]
d <- pivot_longer(d, c("precision", "recall", "F1"), names_to = "metric", values_to = "value")
d$metric <- factor(d$metric, levels = c("precision", "recall", "F1"))
d$variant_type <- factor(d$variant_type, levels = c("SNV", "INDEL"))

p <- ggplot(d, aes(variant_type, value, fill = platform)) +
  geom_col(position = position_dodge(.8), width = .7) +
  geom_text(aes(label = sprintf("%.3f", value)), position = position_dodge(.8),
            vjust = -0.4, size = 2.6) +
  facet_wrap(~metric) +
  scale_fill_manual(values = PLAT) +
  coord_cartesian(ylim = c(0.8, 1.02)) +
  labs(title = "Figure 3 - Cross-platform comparison at nominal 30x (DeepVariant)",
       subtitle = "Same source DNA; mean over three downsampling seeds",
       x = NULL, y = NULL, fill = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.major.x = element_blank())

ggsave(out("Fig3_platform.png"), p, width = 11, height = 4.4, dpi = 300)
cat("->", out("Fig3_platform.png"), "\n")
