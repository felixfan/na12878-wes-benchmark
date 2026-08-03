# Figure 3 — Cross-platform comparison at nominal 30x (DeepVariant). Reads source_data.xlsx.
library(readxl); library(ggplot2); library(dplyr); library(tidyr)
d <- read_excel("../source_data/source_data.xlsx", sheet = "Fig3_platform") |>
  filter(caller == "deepvariant") |>
  pivot_longer(c(precision, recall, F1), names_to = "metric", values_to = "value")
cols <- c(illumina_novaseq6000 = "#4C72B0", bgi_dnbseq_t7 = "#DD8452", genemind_surfseq5000 = "#55A868")
ggplot(d, aes(variant_type, value, fill = platform)) +
  geom_col(position = position_dodge(.8), width = .7) +
  facet_wrap(~factor(metric, c("precision", "recall", "F1"))) +
  coord_cartesian(ylim = c(.8, 1.0)) + scale_fill_manual(values = cols) +
  labs(title = "Figure 3 — Cross-platform comparison at nominal 30x (DeepVariant)", x = NULL, y = NULL) +
  theme_bw() + theme(legend.position = "bottom")
ggsave("Fig3_platform_R.png", width = 12, height = 4.5, dpi = 300)
