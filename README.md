# PTSD proteomic MR analysis code

Current release: **2026.10.20**. This package contains the analysis sources, derived small outputs and report code used for the plasma-protein PTSD study. The default SEM snapshot is the **1,162-protein** result set. Loading the entry point does not run an association analysis.

## Use the supplied results

Extract the complete archive. Keep the folders beside `00_START.R` together. Use the extracted directory on your computer:

```r
source("/path/to/ptsd-proteome-mr-current/00_START.R", encoding="UTF-8")
ptsd_validate()
ptsd_summary()
ptsd_reports(figure2=TRUE, figure4=TRUE, suppfig2=TRUE)
```

Reports read shipped small snapshots and do not need the original GWAS, LD genotypes or GenomicSEM. Figure 2 uses base R and produces PDF plus 600-dpi PNG/TIFF by default. Optional PDF-only export:

```r
options(ptsd.figure2.formats="pdf")
ptsd_reports(figure2=TRUE, figure4=TRUE, suppfig2=TRUE)
```

`ptsd_reports(figure3=TRUE)` additionally plots the original outcome-specific MR coefficients and requires ggplot2, patchwork and ggrepel. `ptsd_reports(figure4=TRUE)` adds the direction-aware base-R brain/blood SMR heatmap. `ptsd_reports(suppfig2=TRUE)` generates the evidence matrix with reverse-estimator and transcript-direction labels. Formatted Additional files 1–3 are supplied separately; reports export component CSVs and plots rather than rewriting Word or Excel submission documents.

## Current SEM results

| Metric | Primary NTOT specification | Source-NEFF sensitivity |
|---|---:|---:|
| Testable proteins | 1,162 | 1,162 |
| Unique instrument SNPs | 1,152 | 1,152 |
| Factor-MR BH-FDR < 0.05 | 51 | 50 |
| Nominal Q_SNP P < 0.05, protein tests | 87 | 85 |
| Q_SNP BH-FDR < 0.05 | 0 | 0 |
| Testable PTSD-prioritized proteins | 9 of 10 | 9 of 10 |
| Factor-MR significant PTSD candidates | 8 of 9 tested | 8 of 9 tested |

UBE2L6, FURIN and CD101 are directly testable. Only SIRPA remains unavailable because rs6136377 is absent from raw MDD GWAS. CD101 has nominal Q_SNP P = 0.00592, adjusted P = 0.563 and factor-MR FDR = 0.236. It provides exploratory heterogeneity evidence, not established PTSD specificity. HYAL1 is the only primary/sensitivity MR-FDR threshold crossing (0.0492 versus 0.0516); all 50 sensitivity hits are primary hits with concordant directions.

The SNP-level analysis uses the genome-wide LDSC trait model. Primary PTSD loading is 0.9261 and standardized residual variance 0.14234. HapMap3 restriction applies to LDSC preparation; the SNP-level analysis uses the expanded, allele/frequency-vetted reference. `results_snapshot/sem_expansion_20261019/` contains the original successful run, full primary/NEFF tables, harmonization audits and control checks. Its `previous_*` files and the older `followup_20261017/sem/` directory are explicitly historical and are not default report inputs.

## Configure an optional analysis rerun

All executable current paths derive from `00_paths.R`. The default project-data directory is the parent of the extracted code folder, so source data can be placed beside the code directory. Reference and protein-data defaults use sibling directories. For another machine, set overrides **before** loading the entry point:

```r
options(ptsd.project_root="/data/PTSD",
        ptsd.path.overrides=list(plink2="/tools/plink2",
             ld_prefix="/data/reference/EUR",
             protein_raw_base="/data/protein_GWAS"))
source("/path/to/code/00_START.R", encoding="UTF-8")
ptsd_check()
```

Use the appropriate GCTA executable for your operating system. `ptsd_check()` checks path existence only. Overrides do not edit release files. `ptsd_validate()` compares issued hashes and parses R files; editing a shipped file deliberately changes its hash.

| Command | Scope |
|---|---|
| `ptsd_run("sem_expand")` | Rebuild expanded SNP statistics and primary/source-NEFF factor MR using audited existing LDSC caches |
| `ptsd_run("sem")` | Rebuild the original primary LDSC model and baseline instrument results |
| `ptsd_run("sem_neff")` | Rebuild the source-NEFF LDSC model and baseline instrument results |
| `ptsd_run("cojo")` | PTSD-only conditional SNP associations |
| `ptsd_run("crp")`, `ptsd_run("coloc")` | cis-CRP MR or uniform dense-region colocalization for all ten proteins |
| `ptsd_run("mediation")`, `ptsd_run("network")` | Exploratory CRP indirect effects or six protein pairs |
| `ptsd_source_audit()`, `ptsd_neff_audit()` | Optional full-file source diagnostics |
| `ptsd_neff_tables()` | Reconstruct S24c from shipped small audit files using base R |
| `ptsd_status()`, `ptsd_originals()` | List current commands or retained supporting sources |

For a full SEM rebuild, first create both baseline LDSC caches with `sem` and `sem_neff`, then run `sem_expand`. The expansion checks compatible parameters, pinned original results, trait effects/SEs, MAF and allele orientation, then verifies eight control fits per specification before fitting additional SNPs. The required LDSC RDS covariance caches are external analysis inputs, not silently fabricated from result tables. If multiple different caches exist, select reviewed directories with `options(ptsd.sem.expand.primary_dir="...", ptsd.sem.expand.neff_dir="...")` or the corresponding configuration fields. Resume is explicit only: `options(ptsd.sem.expand.resume_dir="...")`; the default is fresh extraction. Existing snapshots are never automatically replaced by a new run.

Analyses require their original installed packages (including data.table, R.utils, GenomicSEM, coloc and susieR) and configured source data/executables. No packages are installed automatically. Analyses used Windows R 4.5.2 and GenomicSEM 0.0.5; exact environments are in the run logs.

## Sensitivity annotations

Figure 4 uses signed SMR association strength, displays transcript coefficients and plasma-MR signs, and flags the five qualified opposite-direction tissue/gene results across KHK, SNX18, UBE2L6 and FURIN. Statistical qualification is distinct from direction concordance.

Reverse MR retains six IVW P < 0.05 results. Only CD40 is also nominally significant with the weighted median; CGREF1, FES, KHK, SIRPA and UBE2L6 have IVW-only results. All six are cautionary in Figure S2, and the other four are inconclusive. Forward multi-instrument IVW was tested for SIRPA and CD101; their results are classified as tested and non-significant. Tables S1, S8 and S14 and their notes are shipped as current component CSVs. `checks/check_sensitivity_labels.R` verifies classifications against the displayed source statistics. The historical S2 plotting module is disabled with a pointer to the current renderer.

## Other study results and limits

The main PTSD screen includes 1,258 matched proteins with ten FDR-prioritized candidates; 1,162 refers to common-factor testability. The primary model uses per-SNP PTSD NTOT and sample prevalence 0.1121416; the source-NEFF sensitivity uses provider NEFF and prevalence 0.5. MDD/anxiety settings are fixed. Source diagnostics support the NEFF information scale and a cohort-additive structure, while complete absolute cohort weights and the unique global participant count remain unresolved (S24c). Similar SEM results alone do not validate sample-size definitions, equal effects across disorders or independent replication.

PTSD-only COJO uses 503 EUR participants and 4,686 variants. Conditional P values are 0.0013525 for the FES sentinel and 0.829647 for the FURIN sentinel. These are conditional SNP associations, not conditional protein MR. The FES-tagged colocalization result is prior-sensitive; FURIN lacks strong local PTSD colocalization. No protein pair passes BH-FDR. CRP product estimates do not identify natural indirect effects; FURIN remains unavailable for that analysis because CRP allele frequency is absent.

The affiliation remains **School of Clinical Medicine, Chongqing Medical and Pharmaceutical College, Chongqing, China**.

## Package contents and upload

Upload the extracted contents, preserving this structure, to the GitHub repository and create a release for Zenodo archiving. This ZIP is the complete current code package, not an additional patch to merge manually. `scripts/` contains current analysis and audit modules; `reports/` contains plot/table exports; `results_snapshot/` holds reviewed derived outputs and execution evidence. `retained_modules/` preserves supporting sources for other manuscript analyses; they are optional to execute, have configurable paths/compatibility aliases, and have not all been rerun end-to-end. Keep them when depositing the complete study code; see `OLD_TO_CURRENT.csv` for their scope.

Full raw GWAS/pQTL files, LD genotypes, PLINK/GCTA binaries, credentials and installed packages are excluded. Original logs and manifests retain the execution paths for provenance. No GitHub upload, release publication or Zenodo deposit is performed by the report functions. Use the confirmed version-specific DOI when citing a published deposit; no new DOI is asserted here. `CITATION.cff`, `MANUSCRIPT_MAP.csv`, `SCRIPT_STATUS.csv`, `MD5SUMS.csv` and `SHA256SUMS.csv` describe this release.


## Dense regional colocalization

Dense regional colocalization uses the uniform sentinel-centered GRCh37 +/-500-kb procedure across all ten prioritized proteins. Final coverage is 1,274–3,209 aligned reference-covered SNPs. All ten use p1=p2=1e-4, reference p12=1e-5 and the same seven-value prior grid. No association-P filtering or clumping is performed.

At the reference prior, maximum SuSiE signal-pair H4 is 0.995837 (FES), 0.772537 (CGREF1), 0.471244 (KHK), 0.009622 (FURIN) and 0.002045 (UBE2L6). The maximum is descriptive and need not be the sentinel-containing pair. CGREF1's sentinel is absent from retained protein credible sets. UBE2L6's sentinel-containing pair has H4=1.118e-5. All 462 evaluable signal-pair/prior/N-setting rows are retained.

AKT3, CD40, SNX18 and CD101 have no PTSD 95% credible set under either N setting. CD101 protein fitting is explicitly not run after its two outcome fits converge without a set. SIRPA protein fitting did not converge at 1,000 iterations and is not interpreted. CD40 has high single-causal-variant ABF H4 (0.929543), so its evidence is method-dependent, not a SuSiE-positive finding. Missing credible sets and fitting errors are not negative biological evidence.

`ptsd_reports()` exports S11a–S11d and generates the ten-panel Figure S9 from frozen data. `ptsd_reports(suppfig2=TRUE)` also redraws the current evidence matrix. No analysis rerun is needed to use these results. `ptsd_run("coloc")` is an optional fresh ten-region analysis. Existing results are not overwritten. The implementation evaluates CD101 outcome fits first and skips the protein fit when both outcome fits converge without a credible set. Run-specific execution records and fit statuses are retained with the snapshots.

Multi-instrument MR and MVMR use candidate-window credible-set lead variants; dense-region colocalization is reported separately. The root-level FES/FURIN colocalization snapshots are historical reproducibility controls; current results are in `results_snapshot/coloc_full_region_20261020`. Superseded S11 tables are under `provenance/previous_colocalization_tables`. Verbatim executed scripts are retained as non-entry-point `.R.txt` provenance and contain original local paths; active configuration is relative/overrideable. Public copied logs use path tokens, with source-archive and script hashes retained.
