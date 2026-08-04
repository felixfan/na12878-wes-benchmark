# Figure 6 - CPU versus GPU acceleration of variant calling (132 matched pairs)
# Reads source_data/Fig6_speedup.csv and writes figures/Fig6_speedup.png (300 dpi).
#   Rscript figure_code/Fig6_speedup.R
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

d <- read.csv(src("Fig6_speedup.csv"))
CALLER <- c(gatk4_hc = "#4C72B0", deepvariant = "#DD8452")
LAB <- c(gatk4_hc = "HaplotypeCaller", deepvariant = "DeepVariant")

guides <- data.frame(f = c(5, 10, 20))
p1 <- ggplot(d, aes(call_sec_cpu, call_sec_gpu, colour = caller_pair)) +
  geom_abline(data = guides, aes(slope = 1 / f, intercept = 0),
              linetype = "dashed", colour = "grey60", linewidth = .35) +
  geom_point(size = 1.8, alpha = .85) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = CALLER, labels = LAB) +
  labs(title = "Per-pair calling wall-time", x = "CPU calling time (s)",
       y = "GPU calling time (s)", colour = NULL,
       caption = "dashed guides: 5x, 10x, 20x") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")

lv <- c("5", "10", "15", "30", "50", "100", "200", "full")
med <- d %>%
  mutate(target_depth_x = factor(as.character(target_depth_x), levels = lv)) %>%
  group_by(caller_pair, target_depth_x) %>%
  summarise(speedup = median(speedup), .groups = "drop")

p2 <- ggplot(med, aes(target_depth_x, speedup, colour = caller_pair, group = caller_pair)) +
  geom_line(linewidth = .8) + geom_point(size = 2.2) +
  scale_colour_manual(values = CALLER, labels = LAB) +
  labs(title = "Median speed-up versus depth", x = "Nominal on-target depth",
       y = "Speed-up (CPU / GPU)", colour = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.major.x = element_blank())

p <- (p1 | p2) + plot_annotation(
  title = "Figure 6 - GPU acceleration of variant calling (NVIDIA Parabricks)",
  theme = theme(plot.title = element_text(size = 13, face = "bold")))

ggsave(out("Fig6_speedup.png"), p, width = 11.5, height = 5, dpi = 300)
cat("->", out("Fig6_speedup.png"), "\n")
