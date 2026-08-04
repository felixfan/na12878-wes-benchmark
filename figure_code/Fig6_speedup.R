# Figure 6 - GPU acceleration of variant calling, 132 matched pairs (GRCh38)
# Reads source_data/Fig6_speedup.csv and writes figures/fig6.png (300 dpi).
#   Rscript figure_code/Fig6_speedup.R
# Panels: (a) per-pair wall-time, (b) median speed-up versus depth. No title is drawn on the figure.
# The pairing grid covers 5, 10, 15, 30, 50, 100, 200x and full depth: 300x and 400x were run on CPU
# only (T1 high-depth extension), so they have no GPU twin and cannot appear here.
# Runs from any working directory. Needs: ggplot2, dplyr, patchwork
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(patchwork)})

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

d <- read.csv(src("Fig6_speedup.csv"))
d$caller_pair <- factor(CALLER_LAB[d$caller_pair], levels = names(CALLER_COL))

guides_df <- data.frame(f = c(5, 10, 20))
pa <- ggplot(d, aes(call_sec_cpu, call_sec_gpu, colour = caller_pair)) +
  geom_abline(data = guides_df, aes(slope = 1 / f, intercept = 0),
              linetype = "dashed", colour = "grey60", linewidth = .35) +
  geom_point(size = 1.8, alpha = .85) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = CALLER_COL) +
  labs(x = "CPU calling time (s)", y = "GPU calling time (s)", colour = NULL) +
  theme_bw(base_size = 11)

lv <- c("5", "10", "15", "30", "50", "100", "200", "full")
med <- d %>%
  mutate(target_depth_x = factor(as.character(target_depth_x), levels = lv)) %>%
  group_by(caller_pair, target_depth_x) %>%
  summarise(speedup = median(speedup), .groups = "drop")

pb <- ggplot(med, aes(target_depth_x, speedup, colour = caller_pair, group = caller_pair)) +
  geom_line(linewidth = .8) + geom_point(size = 2.2) +
  scale_colour_manual(values = CALLER_COL) +
  labs(x = "Nominal on-target depth (X)", y = "Speed-up (CPU / GPU wall-time)", colour = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.x = element_blank())

p <- (pa | pb) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(legend.position = "bottom", plot.tag = element_text(face = "bold", size = 12))

ggsave(out("fig6.png"), p, width = 11, height = 4.8, dpi = 300)
cat("->", out("fig6.png"), "\n")
