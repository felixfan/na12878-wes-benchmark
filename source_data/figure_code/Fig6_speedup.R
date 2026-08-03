# Figure 6 — CPU vs GPU acceleration (calling wall-time, 60 pairs). Reads source_data.xlsx.
library(readxl); library(ggplot2); library(dplyr)
d <- read_excel("../source_data/source_data.xlsx", sheet = "Fig6_speedup")
lab <- c(gatk4_hc = "HaplotypeCaller", deepvariant = "DeepVariant")
col <- c(gatk4_hc = "#4C72B0", deepvariant = "#DD8452")
# Panel A: per-pair CPU vs GPU calling time (log-log) with iso-speedup guides
gA <- ggplot(d, aes(call_sec_cpu, call_sec_gpu, color = caller_pair)) +
  geom_abline(slope = 1 / c(5, 10, 20), intercept = 0, linetype = 2, color = "grey70") +
  geom_point(size = 2, alpha = .85) +
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = col, labels = lab) +
  labs(title = "Figure 6A — CPU vs GPU calling time", x = "CPU calling (s)", y = "GPU calling (s)") +
  theme_bw() + theme(legend.position = "bottom")
ggsave("Fig6_speedup_R.png", gA, width = 6, height = 5, dpi = 300)
# Panel B: median speedup vs nominal depth
dk <- function(x) ifelse(x == "full", 1e9, as.numeric(x))
d2 <- d |> group_by(caller_pair, target_depth_x) |>
  summarise(speedup = median(speedup), .groups = "drop") |>
  mutate(ord = dk(target_depth_x)) |> arrange(ord)
gB <- ggplot(d2, aes(reorder(target_depth_x, ord), speedup, color = caller_pair, group = caller_pair)) +
  geom_line(linewidth = 1) + geom_point(size = 2.5) +
  scale_color_manual(values = col, labels = lab) +
  labs(title = "Figure 6B — Speedup vs depth", x = "nominal depth", y = "speedup (CPU/GPU)") +
  theme_bw() + theme(legend.position = "bottom")
ggsave("Fig6_speedup_depth_R.png", gB, width = 6, height = 5, dpi = 300)
