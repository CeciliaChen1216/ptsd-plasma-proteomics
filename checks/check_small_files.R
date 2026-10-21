# Optional check of already-reviewed small outputs; no external GWAS is read.
# First source 00_START.R, then source this file. Uses base R only.
ptsd_check_small_files <- function() {
  P <- getOption("ptsd.current.paths")
  if(is.null(P)) stop("Source 00_START.R first")
  read <- function(p) read.csv(file.path(P$snapshot_dir,p),check.names=FALSE,stringsAsFactors=FALSE)
  primary <- read("factor_MR.csv")
  sensitivity <- read("sem_expansion_20261019/NEFF/factor_MR_expanded.csv")
  stopifnot(nrow(primary)==1162L,nrow(sensitivity)==1162L,
    length(unique(primary$outcome_snp))==1152L,length(unique(sensitivity$outcome_snp))==1152L,
    sum(primary$mr_fdr<.05)==51L,sum(sensitivity$mr_fdr<.05)==50L,
    sum(primary$Q_pval<.05)==87L,sum(sensitivity$Q_pval<.05)==85L,
    sum(primary$Q_fdr<.05)==0L,sum(sensitivity$Q_fdr<.05)==0L)
  candidates <- c("AKT3","CD40","CGREF1","FES","FURIN","SIRPA","CD101","KHK","SNX18","UBE2L6")
  pc <- primary[primary$protein%in%candidates,]
  nc <- sensitivity[sensitivity$protein%in%candidates,]
  stopifnot(nrow(pc)==9L,nrow(nc)==9L,sum(pc$mr_fdr<.05)==8L,sum(nc$mr_fdr<.05)==8L,
    !"SIRPA"%in%primary$protein,primary$mr_fdr[primary$protein=="CD101"]>.05,
    setequal(pc$protein[pc$mr_fdr<.05],nc$protein[nc$mr_fdr<.05]))
  for(d in list(primary,sensitivity)) {
    stopifnot(!anyDuplicated(d$protein),all(is.finite(d$mr_beta)),all(d$mr_se>0),
      max(abs(d$mr_fdr-p.adjust(d$mr_p,"BH")))<1e-12,
      max(abs(d$Q_fdr-p.adjust(d$Q_pval,"BH")))<1e-12,
      max(abs(d$mr_beta-d$beta_factor/d$pqtl_beta))<1e-12,
      max(abs(d$mr_se-d$se_factor/abs(d$pqtl_beta)))<1e-12,
      max(abs(d$mr_p-2*pnorm(-abs(d$mr_beta/d$mr_se))))<1e-12,
      max(abs(d$Q_pval-pchisq(d$Q,df=2,lower.tail=FALSE)))<1e-12)
  }
  stopifnot(setequal(setdiff(primary$protein[primary$mr_fdr<.05],sensitivity$protein[sensitivity$mr_fdr<.05]),"HYAL1"),
    all(sensitivity$protein[sensitivity$mr_fdr<.05]%in%primary$protein[primary$mr_fdr<.05]))
  old <- read("sem_expansion_20261019/primary/previous_stage12_factorMR_results.csv")
  ix <- match(old$protein,primary$protein)
  stopifnot(nrow(old)==391L,!anyNA(ix),max(abs(old$mr_beta-primary$mr_beta[ix]))<1e-12,
    max(abs(old$Q_pval-primary$Q_pval[ix]))<1e-12)
  sirpa <- read("sem_expansion_20261019/candidate_filter_audit.csv")
  stopifnot(!sirpa$MDD_raw_present[sirpa$protein=="SIRPA"])
  cojo <- read("followup_20261017/cojo/COJO_PTSD__COJO_PTSD_comparison.csv")
  stopifnot(nrow(cojo)==2L,all(cojo$status=="estimated"),
    cojo$p_conditional[cojo$target_gene=="FES"]==0.0013525,
    cojo$p_conditional[cojo$target_gene=="FURIN"]==0.829647)
  env <- new.env(parent=globalenv())
  sys.source(file.path(P$package_root,"scripts/PTSD_NEFF_count_patterns.R"),envir=env)
  out <- tempfile("S24c_check_")
  on.exit(unlink(out,recursive=TRUE),add=TRUE)
  x <- env$ptsd_build_neff_tables(P$source_audit_snapshot,out)
  stopifnot(nrow(x$cohorts)==63L,nrow(x$matches)==67L,
    all(x$matches$observed_minus_calculated==0),
    x$diagnostics$count_matched_SNP_rows==2270303,
    x$diagnostics$named_cohorts_in_matches==19L)
  frozen <- read("current_manuscript_components/TableS24c_count_pattern_matches.csv")
  stopifnot(isTRUE(all.equal(x$matches,frozen,check.attributes=FALSE,tolerance=1e-12)))
  cat("PASS: expanded primary/NEFF counts, full-family BH, Wald/Q arithmetic, 391 baseline coefficients, candidate coverage, COJO and reproducible S24c. No association analysis run.\n")
  invisible(TRUE)
}
ptsd_check_small_files()
