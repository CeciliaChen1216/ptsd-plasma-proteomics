# Portable configuration. Loading only records paths; no analysis is run.
# Optional BEFORE 00_START.R:
# options(ptsd.project_root="/data/PTSD", ptsd.path.overrides=list(plink2="/tools/plink2", ld_prefix="/data/EUR"))
package_root <- getOption("ptsd.current.root")
if(is.null(package_root)) stop("Source 00_START.R first")
project_root <- getOption("ptsd.project_root",dirname(package_root))
workspace_root <- dirname(project_root)
P <- list(package_root=package_root, project_root=project_root, workspace_root=workspace_root,
  protein_raw_base=file.path(workspace_root,"Longcovid/data/raw/pqtl_ukb_ppp/cis_regions"),
  plink2=file.path(workspace_root,"Longcovid/tools",if(.Platform$OS.type=="windows") "plink2.exe" else "plink2"),
  gcta=file.path(project_root,"tools/GCTA_1.95.1/gcta-1.95.1-windows-x86_64/bin/gcta64.exe"),
  ld_prefix=file.path(workspace_root,"reference/1000G_v3/EUR"),
  hm3=file.path(workspace_root,"reference/ldsc/w_hm3.snplist"),
  ldsc=file.path(workspace_root,"reference/ldsc/eur_w_ld_chr"),
  sumstats_reference=file.path(workspace_root,"reference/reference.custom.maf.0.005.txt.gz"))
P$result_dir <- file.path(project_root,"results")
P$pqtl_rds <- file.path(project_root,"pQTL/ukbppp_cis_pqtl.rds")
P$ptsd_ma <- file.path(P$result_dir,"ptsd_freeze3_smr.ma")
P$main_mr_rds <- file.path(P$result_dir,"mr_all_outcomes.rds")
P$crp_gwas <- file.path(project_root,"validation_data/external_data/crp_gwas/GCST90029070_buildGRCh37.tsv.gz")
P$source_input_dir <- file.path(project_root,"data_source_check_20261016")
P$source_casecontrol <- file.path(P$source_input_dir,"eur_ptsdcasecontrol_pcs_v4_aug3_2021.vcf.gz")
P$mdd_gwas <- file.path(project_root,"GWAS/PGC_MDD_2025/27061255/daner/daner_pgc_mdd_no23andMe-noUKBB_eur_hg19_v3.49.24.11.neff.gz")
P$anx_gwas <- file.path(project_root,"GWAS/PGC_Anxiety_2026/31389910/ANX_2026_daner_fullANX_v12_woUTAH_11022026.gz")
P$snapshot_dir <- file.path(package_root,"results_snapshot")
P$coloc_snapshot <- file.path(P$snapshot_dir,"coloc_full_region_20261020")
P$followup_snapshot <- file.path(P$snapshot_dir,"followup_20261017")
P$sem_expansion_snapshot <- file.path(P$snapshot_dir,"sem_expansion_20261019")
P$sem_baseline_snapshot <- file.path(P$sem_expansion_snapshot,"primary/previous_stage12_factorMR_results.csv")
P$crp_step2_results <- file.path(P$snapshot_dir,"CRP_to_PTSD_results.csv")
P$ptsd_region <- file.path(P$snapshot_dir,"chr15_PTSD_with_sample_sizes.csv")
P$source_audit_snapshot <- file.path(P$snapshot_dir,"source_audit_20261017")
P$report_base <- file.path(project_root,"reports_current_20261020")
P$report_dir <- P$report_base
# Optional cache overrides; the expansion otherwise audits compatible caches under result_dir.
P$sem_primary_cache <- ""
P$sem_neff_cache <- ""
overrides <- getOption("ptsd.path.overrides",list())
if(!is.list(overrides)||is.null(names(overrides))&&length(overrides)) stop("ptsd.path.overrides must be a named list")
unknown <- setdiff(names(overrides),names(P))
if(length(unknown)) stop("Unknown path override: ",paste(unknown,collapse=", "))
P[names(overrides)] <- overrides
options(ptsd.current.paths=P,ptsd.project_root=project_root)
