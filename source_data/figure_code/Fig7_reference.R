# Figure 7 — Reference-build comparison (b37 vs GRCh38, full depth). Reads source_data.xlsx.
library(readxl); library(ggplot2); library(dplyr)
d <- read_excel("../source_data/source_data.xlsx", sheet = "Fig7_reference") |>
  mutate(platform = sub("_.*", "", platform))
ggplot(d, aes(platform, F1, fill = reference)) +
  geom_col(position = position_dodge(.8), width = .7) +
  geom_text(aes(label = sprintf("%.3f", F1)), position = position_dodge(.8), vjust = -.3, size = 2.6) +
  facet_grid(variant_type ~ caller) +
  coord_cartesian(ylim = c(.75, 1.02)) +
  scale_fill_manual(values = c(b37 = "#8C8C8C", grch38 = "#4C72B0")) +
  labs(title = "Figure 7 — Reference-build comparison (full depth)", x = NULL, y = "F1") +
  theme_bw() + theme(legend.position = "bottom")
ggsave("Fig7_reference_R.png", width = 10.5, height = 7.8, dpi = 300)
