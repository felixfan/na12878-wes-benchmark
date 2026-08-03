# Figure 5 — CPU vs GPU F1 concordance (60 matched pairs). Reads source_data.xlsx.
library(readxl); library(ggplot2)
d <- read_excel("../source_data/source_data.xlsx", sheet = "Fig5_concordance")
lo <- min(d$CPU_F1, d$GPU_F1) - .004
ggplot(d, aes(CPU_F1, GPU_F1, color = caller_pair)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
  geom_point(size = 2.4, alpha = .85) +
  coord_equal(xlim = c(lo, 1.001), ylim = c(lo, 1.001)) +
  labs(title = "Figure 5 — CPU vs GPU concordance (F1)", x = "CPU F1", y = "GPU F1") +
  theme_bw() + theme(legend.position = c(.28, .9))
ggsave("Fig5_concordance_R.png", width = 5.6, height = 5.4, dpi = 300)
