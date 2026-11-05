# Release history

## Version 1.1.0

Prepared on 2026-11-05.

- Replaced the 66 formal GPU HaplotypeCaller runs with interval-restricted BQSR results. Both Parabricks BQSR steps use the V6 BED and zero additional interval padding.
- Replaced 792 benchmark rows and the corresponding run manifest and timing records. The archive contains 5,448 rows and 454 run identifiers, with 450 formal runs.
- Recomputed the 66 HaplotypeCaller site-concordance and quality comparisons. The archive retains 134 pairs, of which 132 are used in the analyses.
- Updated Figure 5, Figure 6, Table3, Table4 and Supplementary Table S1 source data. Regenerated the figure PNG and PDF files.
- Added per-run protocol, target BED checksum, stage timings and previous affected records under results/revisions/v1.1.0.
- Added protocol-aware GPU BQSR caches and a guard against using incompatible existing run outputs.
- Added a SHA-256 file inventory.
- Added Git attributes that preserve the bytes used by the file inventory across Windows and Linux checkouts.

The four acceptance runs and their two matched pairs remain historical, excluded records. Their original processing provenance is preserved in the revision data.

The HaplotypeCaller median calling-time speedup is approximately 9.9. Its SNV and INDEL genotype concordance percentages are 99.999 and 99.993. Discordant-site quality medians are based on 91 nonempty pairs across both callers.

## Version 1.0.1

The previous published snapshot is archived at [Zenodo](https://doi.org/10.5281/zenodo.21792675). Its formal GPU HaplotypeCaller BQSR steps used the input WES BAM without an interval restriction. That configuration is superseded by version 1.1.0 for formal comparisons.
