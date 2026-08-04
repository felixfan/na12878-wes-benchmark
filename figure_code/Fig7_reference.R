# Figure 7 - Reference-build comparison, GRCh38 versus GRCh37/b37 (full depth)
# Reads source_data/Fig7_reference.csv and writes figures/Fig7_reference.png (300 dpi).
#   Rscript figure_code/Fig7_reference.R
# Runs from any working directory. Needs: ggplot2
suppressPackageStartupMessages(library(ggplot2))

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

d <- read.csv(src("Fig7_reference.csv"))
d$variant_type <- factor(d$variant_type, levels = c("SNV", "INDEL"))
d$reference <- factor(d$reference, levels = c("b37", "grch38"))
d$platform_short <- sub("_.*", "", d$platform)

p <- ggplot(d, aes(platform_short, F1, fill = reference)) +
  geom_col(position = position_dodge(.75), width = .65) +
  geom_text(aes(label = sprintf("%.3f", F1)), position = position_dodge(.75),
            vjust = -0.4, size = 2.5) +
  facet_grid(variant_type ~ caller) +
  scale_fill_manual(values = c(b37 = "#8C8C8C", grch38 = "#4C72B0"),
                    labels = c(b37 = "GRCh37/b37", grch38 = "GRCh38")) +
  coord_cartesian(ylim = c(0.75, 1.04)) +
  labs(title = "Figure 7 - Reference-build comparison at full depth",
       subtitle = "Same reads and callers, evaluated against the GIAB truth set of the matching build",
       x = NULL, y = "F1", fill = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.major.x = element_blank())

ggsave(out("Fig7_reference.png"), p, width = 9.5, height = 7, dpi = 300)
cat("->", out("Fig7_reference.png"), "\n")
