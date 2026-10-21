# One entry point. Loading defines functions only; it does not run analyses.
.ptsd_sources <- vapply(sys.frames(), function(e) {
  f <- e$ofile
  if(is.character(f) && length(f)==1L) f else NA_character_
}, character(1))
.ptsd_sources <- .ptsd_sources[!is.na(.ptsd_sources)]
.ptsd_sources <- .ptsd_sources[basename(.ptsd_sources)=="00_START.R"]
if(!length(.ptsd_sources)) stop("Use source('.../00_START.R', encoding='UTF-8')")
.ptsd_file <- tail(.ptsd_sources,1L)
.ptsd_root <- dirname(normalizePath(.ptsd_file,winslash="/",mustWork=TRUE))
options(ptsd.current.root=.ptsd_root)
source(file.path(.ptsd_root,"00_paths.R"),local=TRUE,encoding="UTF-8")

ptsd_status <- function(all=FALSE) {
  P <- getOption("ptsd.current.paths")
  d <- read.csv(file.path(P$package_root,"SCRIPT_STATUS.csv"),stringsAsFactors=FALSE)
  if(!all) d <- d[d$scope %in% c("current","audit","report"),]
  print(d,row.names=FALSE)
  invisible(d)
}

ptsd_check <- function() {
  P <- getOption("ptsd.current.paths")
  paths <- c(pqtl=P$pqtl_rds,main_MR=P$main_mr_rds,CRP=P$crp_gwas,PTSD_ma=P$ptsd_ma,
    source_casecontrol=P$source_casecontrol,MDD=P$mdd_gwas,ANX=P$anx_gwas,
    PLINK2=P$plink2,GCTA=P$gcta,LD_bed=paste0(P$ld_prefix,".bed"),
    LD_bim=paste0(P$ld_prefix,".bim"),LD_fam=paste0(P$ld_prefix,".fam"),
    HM3=P$hm3,sumstats_reference=P$sumstats_reference,
    protein_directory=P$protein_raw_base,frozen_CRP_result=P$crp_step2_results,
    PTSD_region=P$ptsd_region)
  d <- data.frame(item=names(paths),exists=file.exists(paths),path=unname(paths),row.names=NULL)
  print(d,row.names=FALSE)
  cat("Path check only. GCTA is needed for cojo; SEM reference files are needed only for SEM. No analysis run.\n")
  invisible(d)
}

ptsd_validate <- function() {
  P <- getOption("ptsd.current.paths")
  m <- read.csv(file.path(P$package_root,"MD5SUMS.csv"),stringsAsFactors=FALSE)
  files <- file.path(P$package_root,m$relative_path)
  observed <- unname(tools::md5sum(files))
  d <- data.frame(file=m$relative_path,exists=file.exists(files),
    matches=!is.na(observed)&observed==m$md5,stringsAsFactors=FALSE)
  print(d[!d$exists | !d$matches,,drop=FALSE],row.names=FALSE)
  if(any(!d$exists | !d$matches)) stop("Bundle differs from issued files. If you edited paths, that file should differ.")
  current_R <- c(file.path(P$package_root,c("00_START.R","00_paths.R","00_config.R")),
    list.files(file.path(P$package_root,"scripts"),pattern="[.]R$",full.names=TRUE),
    list.files(file.path(P$package_root,"reports"),pattern="[.]R$",full.names=TRUE))
  invisible(lapply(current_R,parse,encoding="UTF-8"))
  cat("All ",nrow(d)," listed files match this release; current R scripts parse. No analysis run.\n",sep="")
  invisible(d)
}

ptsd_run <- function(module) {
  map <- c(sem="PTSD_GenomicSEM_clean_rerun.R",crp="PTSD_CRP_cis_MR_rerun_v2.R",
    coloc="PTSD_coloc_full_region_rerun.R",network="PTSD_network_MR_rerun_v2.R",
    mediation="PTSD_CRP_mediation_rerun_v2.R",sem_expand="PTSD_SEM_expand_instruments.R")
  choices <- c(names(map),"sem_neff","cojo")
  if(length(module)!=1L || !module %in% choices) stop("Choose: ",paste(choices,collapse=", "))
  P <- getOption("ptsd.current.paths")
  old <- getwd(); on.exit(setwd(old),add=TRUE); setwd(P$package_root)
  dir.create(P$result_dir,recursive=TRUE,showWarnings=FALSE)
  env <- new.env(parent=globalenv())
  cat("Requested: ",module,". New timestamped output; reviewed snapshots stay unchanged.\n",sep="")
  if(module %in% c("sem_neff","cojo")) {
    opt <- options(ptsd.followup.package_root=P$package_root,
      ptsd.followup.mode=if(module=="sem_neff") "sem" else "cojo",
      ptsd.followup.load_only=TRUE)
    on.exit(options(opt),add=TRUE)
    if(nzchar(P$gcta)) {
      opt_gcta <- options(ptsd.followup.gcta=P$gcta)
      on.exit(options(opt_gcta),add=TRUE)
    }
    sys.source(file.path(P$package_root,"scripts/PTSD_SEM_COJO_followup_20261016.R"),envir=env)
    return(invisible(env$ptsd_followup_20261016()))
  }
  if(module=="sem") cat("Prepare the recorded primary LDSC/baseline cache; current expanded SNP results require sem_expand afterward.\n")
  source(file.path(P$package_root,"scripts",map[[module]]),local=env,encoding="UTF-8")
}

ptsd_source_audit <- function() {
  P <- getOption("ptsd.current.paths")
  opt <- options(ptsd.source_check.load_only=TRUE); on.exit(options(opt),add=TRUE)
  env <- new.env(parent=globalenv())
  sys.source(file.path(P$package_root,"scripts/PTSD_source_check_20261016.R"),envir=env)
  invisible(env$ptsd_source_check(input_dir=P$source_input_dir,project_dir=P$project_root))
}

ptsd_reports <- function(figure2=FALSE, figure3=FALSE, figure4=FALSE, suppfig2=FALSE) {
  P <- getOption("ptsd.current.paths")
  dir.create(P$report_base,recursive=TRUE,showWarnings=FALSE)
  P$report_dir <- tempfile(paste0("export_",format(Sys.time(),"%Y%m%d_%H%M%S"),"_"),tmpdir=P$report_base)
  if(!dir.create(P$report_dir)) stop("Cannot create report directory")
  opt <- options(ptsd.current.paths=P); on.exit(options(opt),add=TRUE)
  for(n in c("make_reviewed_reports.R","make_followup_reports.R"))
    source(file.path(P$package_root,"reports",n),local=new.env(parent=globalenv()),encoding="UTF-8")
  if(figure2) source(file.path(P$package_root,"reports/make_figure2_current.R"),
    local=new.env(parent=globalenv()),encoding="UTF-8")
  if(figure3) source(file.path(P$package_root,"reports/make_figure3_current.R"),
    local=new.env(parent=globalenv()),encoding="UTF-8")
  if(figure4) source(file.path(P$package_root,"reports/make_figure4_current.R"),
    local=new.env(parent=globalenv()),encoding="UTF-8")
  if(suppfig2) source(file.path(P$package_root,"reports/make_suppfig2_current.R"),
    local=new.env(parent=globalenv()),encoding="UTF-8")
  cat("Reviewed snapshots exported to: ",P$report_dir,"\nNo GWAS analysis run; Word/XLSX manuscripts not rewritten.\n",sep="")
  invisible(P$report_dir)
}

ptsd_neff_audit <- function() {
  P <- getOption("ptsd.current.paths")
  opt <- options(ptsd.neff_source_audit.load_only=TRUE); on.exit(options(opt),add=TRUE)
  env <- new.env(parent=globalenv())
  sys.source(file.path(P$package_root,"scripts/PTSD_NEFF_source_audit_20261017.R"),envir=env)
  invisible(env$ptsd_neff_source_audit(source_file=P$source_casecontrol,
                                      output_parent=P$source_input_dir))
}

ptsd_neff_tables <- function(audit_dir=NULL, output_dir=NULL) {
  P <- getOption("ptsd.current.paths")
  if(is.null(audit_dir)) audit_dir <- P$source_audit_snapshot
  if(is.null(output_dir)) {
    dir.create(P$report_base,recursive=TRUE,showWarnings=FALSE)
    output_dir <- tempfile("neff_tables_",tmpdir=P$report_base)
  }
  env <- new.env(parent=globalenv())
  sys.source(file.path(P$package_root,"scripts/PTSD_NEFF_count_patterns.R"),envir=env)
  result <- env$ptsd_build_neff_tables(audit_dir=audit_dir,output_dir=output_dir)
  cat("Table S24c regenerated from small audit files: ",output_dir,
      "\nNo full-GWAS read or association analysis.\n",sep="")
  invisible(result)
}

ptsd_summary <- function() {
 P <- getOption("ptsd.current.paths")
 x <- read.csv(file.path(P$sem_expansion_snapshot,"summary_counts_all_modes.csv"),stringsAsFactors=FALSE)
 print(x,row.names=FALSE)
 cat("Nine of ten PTSD candidates are testable; eight pass factor-MR FDR in both specifications.\n",
     "SIRPA is absent from raw MDD GWAS. CD101 is testable but not factor-MR FDR significant.\n",sep="")
 invisible(x)
}

ptsd_originals <- function() {
  P <- getOption("ptsd.current.paths")
  n <- list.files(file.path(P$package_root,"retained_modules"),pattern="[.]R$")
  print(n)
  cat("Historical reference modules, not all revalidated. See OLD_TO_CURRENT.csv; do not run the old COJO/table appenders for the updated manuscript.\n")
  invisible(n)
}

cat("PTSD current bundle 2026-10-20 revision 1 loaded. NO analysis run.\n",
  "ptsd_check()               : local input paths\n",
  "ptsd_summary()             : reviewed expanded SEM counts (1162 proteins)\n",
  "ptsd_validate()            : release file integrity\n",
  "ptsd_status()              : current module list\n",
  "ptsd_reports()             : export reviewed tables/plots, no GWAS rerun\n",
  "ptsd_reports(figure2=TRUE)  : also generate recorded-primary Figure 2\n",
  "ptsd_reports(figure3=TRUE)  : Figure 3 from frozen full MR tables\n",
  "ptsd_reports(figure4=TRUE)  : direction-aware brain/blood SMR heatmap\n",
  "ptsd_reports(suppfig2=TRUE) : current evidence matrix including dense colocalization\n",
  "ptsd_neff_tables()         : regenerate S24c from frozen small audit files\n",
  "ptsd_neff_audit()          : OPTIONAL full-file source-NEFF audit\n",
  "ptsd_run('sem_expand')     : OPTIONAL expanded SNP analysis under both reviewed LDSC specifications\n",
  "ptsd_run('sem') / sem_neff  : OPTIONAL baseline LDSC/cache preparation; also crp/coloc/mediation/network/cojo\n",
  "ptsd_source_audit()        : OPTIONAL read-only source audit\n",
  "ptsd_originals()           : list historical reference modules\n",sep="")
