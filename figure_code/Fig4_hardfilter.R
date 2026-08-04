# Figure 4 - GATK INDEL hard-filtering, before and after (GATK4 HaplotypeCaller, full depth, GRCh38)
# Reads source_data/Fig4_hardfilter.csv and writes figures/Fig4_hardfilter.png (300 dpi).
#   Rscript figure_code/Fig4_hardfilter.R
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

d <- read.csv(src("Fig4_hardfilter.csv"), check.names = FALSE)
d <- pivot_longer(d, c("F1", "precision", "recall"), names_to = "metric", values_to = "value")
d$metric <- factor(d$metric, levels = c("F1", "precision", "recall"))
d$variant_type <- factor(d$variant_type, levels = c("SNV", "INDEL"))
d$filter <- factor(d$filter, levels = c("raw (unfiltered)", "hard-filtered"))

p <- ggplot(d, aes(variant_type, value, fill = filter)) +
  geom_col(position = position_dodge(.75), width = .65) +
  geom_text(aes(label = sprintf("%.3f", value)), position = position_dodge(.75),
            vjust = -0.4, size = 2.8) +
  facet_wrap(~metric) +
  scale_fill_manual(values = c("raw (unfiltered)" = "#B0B0B0", "hard-filtered" = "#4C72B0")) +
  coord_cartesian(ylim = c(0.6, 1.05)) +
  labs(title = "Figure 4 - GATK4 HaplotypeCaller: INDEL-only hard-filtering, full depth",
       subtitle = "SNVs are deliberately left unfiltered, so the SNV bars are identical by construction",
       x = NULL, y = NULL, fill = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.major.x = element_blank())

ggsave(out("Fig4_hardfilter.png"), p, width = 10.5, height = 4.4, dpi = 300)
cat("->", out("Fig4_hardfilter.png"), "\n")
