# Export already-reviewed frozen data. No refitting or automatic adoption of
# the newest local analysis; a future rerun needs review before snapshot update.
(function() {
  P <- getOption("ptsd.current.paths")
  if(is.null(P)) stop("Source 00_START.R first")
  files <- c(
    COJO_PTSD_comparison_20261017.csv="followup_20261017/cojo/COJO_PTSD__COJO_PTSD_comparison.csv",
    SEM_source_NEFF_full_results.csv="sem_expansion_20261019/NEFF/factor_MR_expanded.csv",
    SEM_source_NEFF_comparison.csv="sem_expansion_20261019/comparison_primary_NEFF.csv",
    TableS1_current.csv="current_manuscript_components/TableS1.csv",
    TableS8_current.csv="current_manuscript_components/TableS8.csv",
    TableS11a_current.csv="current_manuscript_components/TableS11a.csv",
    TableS11b_current.csv="current_manuscript_components/TableS11b.csv",
    TableS11c_current.csv="current_manuscript_components/TableS11c.csv",
    TableS11d_current.csv="current_manuscript_components/TableS11d.csv",
    TableS14_current.csv="current_manuscript_components/TableS14.csv",
    TableS3a_current.csv="current_manuscript_components/TableS3a.csv",
    TableS3b_current.csv="current_manuscript_components/TableS3b.csv",
    TableS3c_current.csv="current_manuscript_components/TableS3c.csv",
    TableS22_current.csv="current_manuscript_components/TableS22.csv",
    TableS24b_current.csv="current_manuscript_components/TableS24b.csv",
    TableS26b_current.csv="current_manuscript_components/TableS26b.csv",
    TableS25_current.csv="current_manuscript_components/TableS25.csv",
    TableS26_current.csv="current_manuscript_components/TableS26.csv",
    SEM_summary_counts.csv="current_manuscript_components/SEM_summary_counts.csv",
    SEM_candidate_primary.csv="sem_expansion_20261019/primary/candidate_results.csv",
    SEM_candidate_source_NEFF.csv="sem_expansion_20261019/NEFF/candidate_results.csv")
  src <- file.path(P$snapshot_dir,unname(files))
  if(any(!file.exists(src))) stop("Missing reviewed snapshot: ",paste(src[!file.exists(src)],collapse=", "))
  ok <- file.copy(src,file.path(P$report_dir,names(files)),overwrite=FALSE)
  if(!all(ok)) stop("Some reviewed components could not be exported")
  env <- new.env(parent=globalenv())
  sys.source(file.path(P$package_root,"scripts/PTSD_NEFF_count_patterns.R"),envir=env)
  env$ptsd_build_neff_tables(P$source_audit_snapshot,P$report_dir)
  writeLines(c("Frozen reviewed components corresponding to the manuscript and Oct20 dense-colocalization code release.",
    "S24c is recalculated from the small full-file audit outputs, not from a full GWAS read.",
    "S22: conditional SNP-to-PTSD log OR, not protein-level MR.",
    "S24b/S26b: source-NEFF sensitivity, not an automatically replaced primary model.",
    "Reports do not edit the complete supplementary workbook or Word documents."),
    file.path(P$report_dir,"FOLLOWUP_COMPONENTS_README.txt"))
})()
