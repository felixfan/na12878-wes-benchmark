# Figure 2 — Depth-saturation curves (three platforms, DeepVariant). Reads source_data.xlsx.
library(readxl); library(ggplot2); library(dplyr)
d <- read_excel("../source_data/source_data.xlsx", sheet = "Fig2_saturation") |>
  filter(caller == "deepvariant")
cols <- c(illumina_novaseq6000 = "#4C72B0", bgi_dnbseq_t7 = "#DD8452", genemind_surfseq5000 = "#55A868")
ggplot(d, aes(mean_on_target_depth, F1_mean, color = platform)) +
  geom_line(linewidth = .8) + geom_point(size = 1.8) +
  geom_errorbar(aes(ymin = F1_mean - F1_std, ymax = F1_mean + F1_std), width = .03) +
  facet_wrap(~variant_type) +
  scale_x_log10(breaks = c(8, 10, 15, 20, 30, 50, 75, 100, 200, 400, 700)) +
  scale_color_manual(values = cols) +
  labs(title = "Figure 2 — Depth-saturation of variant-calling accuracy (DeepVariant)",
       x = "On-target mean depth (X)", y = "F1") +
  theme_bw() + theme(legend.position = "bottom")
ggsave("Fig2_saturation_R.png", width = 11, height = 4.8, dpi = 300)
