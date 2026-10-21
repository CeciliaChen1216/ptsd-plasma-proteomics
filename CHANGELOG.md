# Release 2026.10.20

This release contains analysis code, derived result snapshots and figure/table exporters for the plasma-protein PTSD study.

- The primary PTSD screen includes 1,258 proteins and ten FDR-prioritized candidates.
- Genomic SEM includes 1,162 protein tests representing 1,152 SNPs. Factor-MR identifies 51 primary and 50 source-NEFF FDR associations. Nine prioritized proteins are testable and eight pass factor-MR FDR correction in both specifications. No QSNP test passes FDR correction.
- Dense regional colocalization uses sentinel-centered GRCh37 +/-500-kb windows and uniform QC and priors for all ten prioritized proteins. All evaluable signal pairs and fit-status distinctions are supplied. Missing credible sets and nonconvergence are not interpreted as negative biological evidence.
- Figure 4 distinguishes transcript-effect direction from plasma-protein MR direction. Figure S2 distinguishes estimator agreement, IVW-only signals, tested non-significance and unavailable results.
- Multi-instrument MR and MVMR use candidate-window credible-set lead variants. The conditional chromosome 15 analyses report SNP associations, not conditional protein MR.
- `checks/audit_primary_coverage_and_MDD.py` reads the supplied snapshots, verifies coverage and factor/MDD comparison counts, and does not write files.
- `00_START.R` loads functions without running an association analysis. `ptsd_validate()` checks issued MD5 values and R syntax; optional analyses require configured external input data and executables.

To export supplied results, load `00_START.R` and run `ptsd_reports(figure2=TRUE, figure4=TRUE, suppfig2=TRUE)`. Full raw GWAS/pQTL data, reference genotypes, executables and installed packages are not redistributed. Scientific assumptions and limitations are described in README.md and the accompanying manuscript methods.
