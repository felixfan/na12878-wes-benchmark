# Figure 4 — GATK INDEL hard-filtering before/after (GATK4, full depth). Reads source_data.xlsx.
library(readxl); library(ggplot2); library(dplyr); library(tidyr)
d <- read_excel("../source_data/source_data.xlsx", sheet = "Fig4_hardfilter") |>
  pivot_longer(c(F1, precision), names_to = "metric", values_to = "value")
ggplot(d, aes(variant_type, value, fill = filter)) +
  geom_col(position = position_dodge(.7), width = .6) +
  geom_text(aes(label = sprintf("%.3f", value)), position = position_dodge(.7), vjust = -.3, size = 3) +
  facet_wrap(~metric) + coord_cartesian(ylim = c(.6, 1.03)) +
  scale_fill_manual(values = c(`raw (unfiltered)` = "#B0B0B0", `hard-filtered` = "#4C72B0")) +
  labs(title = "Figure 4 — GATK INDEL hard-filtering before/after (full depth)", x = NULL, y = NULL) +
  theme_bw() + theme(legend.position = "bottom")
ggsave("Fig4_hardfilter_R.png", width = 10.5, height = 4.8, dpi = 300)
