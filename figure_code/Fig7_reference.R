# Figure 7 - Reference-build comparison, GRCh37/b37 versus GRCh38 (full depth)
# Reads source_data/Fig7_reference.csv and writes figures/fig7.png (300 dpi).
#   Rscript figure_code/Fig7_reference.R
# Panels: (a) SNV, GATK4 HaplotypeCaller  (b) SNV, DeepVariant
#         (c) INDEL, GATK4 HaplotypeCaller (d) INDEL, DeepVariant. No title is drawn on the figure.
# Runs from any working directory. Needs: ggplot2, patchwork
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})

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

PLAT_KEY <- c("illumina_novaseq6000", "bgi_dnbseq_t7", "genemind_surfseq5000")
PLAT_LAB <- c("Illumina NovaSeq 6000", "BGI DNBSEQ-T7", "GeneMind SURFSeq 5000")
PLAT_COL <- setNames(c("#4C72B0", "#DD8452", "#55A868"), PLAT_LAB)
REF_LAB <- c(b37 = "GRCh37/b37", grch38 = "GRCh38")

d <- read.csv(src("Fig7_reference.csv"))
d$platform <- factor(PLAT_LAB[match(d$platform, PLAT_KEY)], levels = PLAT_LAB)
d$reference <- factor(REF_LAB[d$reference], levels = REF_LAB)

panel <- function(vt, caller, ylab) {
  s <- d[d$variant_type == vt & d$caller == caller, ]
  ggplot(s, aes(reference, F1, fill = platform)) +
    geom_col(position = position_dodge(.8), width = .7) +
    geom_text(aes(label = sprintf("%.3f", F1)), position = position_dodge(.8),
              vjust = -0.4, size = 2.4) +
    scale_fill_manual(values = PLAT_COL) +
    coord_cartesian(ylim = c(0.75, 1.04)) +
    labs(x = NULL, y = ylab, fill = NULL) +
    theme_bw(base_size = 11) +
    theme(panel.grid.major.x = element_blank())
}

p <- (panel("SNV", "gatk4_hc", "SNV F1") | panel("SNV", "deepvariant", "SNV F1")) /
     (panel("INDEL", "gatk4_hc", "INDEL F1") | panel("INDEL", "deepvariant", "INDEL F1")) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(legend.position = "bottom", plot.tag = element_text(face = "bold", size = 12))

ggsave(out("fig7.png"), p, width = 9.5, height = 7, dpi = 300)
cat("->", out("fig7.png"), "\n")
