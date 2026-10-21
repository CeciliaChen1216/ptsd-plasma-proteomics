P <- getOption("ptsd.current.paths")
if(is.null(P)) stop("Source 00_START.R first")
# PTSD source-NEFF diagnostic, 2026-10-17.
# Read-only audit of the verified official case/control GWAS.
# Run this file with source(..., encoding="UTF-8").
# It creates one NEW audit directory; it does not run MR/LDSC/SEM/COJO,
# alter the current code bundle, or replace any study result.
# The SE-based effective N is an approximation, not proof of the source NEFF definition.

ptsd_neff_source_audit <- function(
  source_file=file.path(P$project_root,"data_source_check_20261016/eur_ptsdcasecontrol_pcs_v4_aug3_2021.vcf.gz"),
  output_parent=file.path(P$project_root,"data_source_check_20261016")
) {
  needed_packages <- c("data.table","R.utils")
  absent <- needed_packages[!vapply(needed_packages,requireNamespace,logical(1),quietly=TRUE)]
  if(length(absent)) stop("Missing installed packages: ",paste(absent,collapse=", "))
  suppressPackageStartupMessages(library(data.table))
  if(!file.exists(source_file)) stop("Source file not found: ",source_file)
  if(!dir.exists(output_parent)) stop("Output parent not found: ",output_parent)
  out <- tempfile(paste0("neff_source_audit_",format(Sys.time(),"%Y%m%d_%H%M%S"),"_"),tmpdir=output_parent)
  if(!dir.create(out)) stop("Cannot create a new audit directory")
  sendback <- file.path(out,"SEND_BACK");dir.create(sendback)
  old_sink <- sink.number();old_message_sink <- sink.number(type="message")
  con <- file(file.path(out,"RUN.log"),"wt")
  sink(con,split=TRUE);sink(con,type="message")
  on.exit({
    if(sink.number(type="message")!=old_message_sink) sink(type="message")
    while(sink.number()>old_sink) sink()
    close(con)
    small <- list.files(out,pattern="[.](csv|txt|log)$",full.names=TRUE)
    if(length(small)) file.copy(small,sendback,overwrite=FALSE)
    cat("\nAudit output: ",out,"\nSend the small files in: ",sendback,"\n",sep="")
  },add=TRUE)
  tryCatch({
    cat("STEP 1/3: Verify the official source\n")
    expected_md5 <- "bf0fb5b16bf91e640af14784e9d3a829"
    md5 <- unname(tools::md5sum(source_file))
    write.csv(data.frame(source_file=source_file,bytes=file.info(source_file)$size,
      observed_md5=md5,expected_md5=expected_md5),file.path(out,"source_manifest.csv"),row.names=FALSE)
    if(!identical(md5,expected_md5)) stop("MD5 differs from the verified Figshare v3 case/control release; inspect before proceeding")
    gz <- gzfile(source_file,"rt");header <- character()
    tryCatch({repeat {
      z <- readLines(gz,n=1,warn=FALSE)
      if(!length(z)) stop("No #CHROM header")
      header <- c(header,z)
      if(startsWith(z,"#CHROM")) break
      if(length(header)>10000L) stop("Unexpected source header")
    }},finally=close(gz))
    writeLines(header,file.path(out,"source_header.txt"))
    writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
    fields <- c("#CHROM","ID","A1","A2","FREQ","BETA","SE","NCASE","NCON","NEFF","NTOT")
    d <- fread(source_file,skip="#CHROM",select=fields,showProgress=TRUE)
    setnames(d,"#CHROM","CHR")
    if(nrow(d)!=6545577L) stop("Source row count differs from the audited release")
    if(anyDuplicated(d$ID)) stop("Duplicate SNP IDs: inspect source before further interpretation")
    cat("Source checksum matched; rows: ",nrow(d),"\n",sep="")

    cat("STEP 2/3: Counts and SE-based diagnostic (no statistical model fitting)\n")
    summarize <- function(x,name,subset) {
      x <- as.numeric(x);x <- x[is.finite(x)]
      if(!length(x)) return(data.table(subset=subset,quantity=name,n=0L,
        min=NA_real_,q01=NA_real_,q05=NA_real_,median=NA_real_,q95=NA_real_,q99=NA_real_,max=NA_real_))
      q <- quantile(x,c(.01,.05,.5,.95,.99),names=FALSE)
      data.table(subset=subset,quantity=name,n=length(x),min=min(x),q01=q[1],q05=q[2],
        median=q[3],q95=q[4],q99=q[5],max=max(x))
    }
    field_summary <- rbindlist(lapply(c("NCASE","NCON","NEFF","NTOT"),
      function(nm) summarize(d[[nm]],nm,"all_source_rows")))
    fwrite(field_summary,file.path(out,"source_sample_size_summary.csv"))
    checks <- data.table(check=c("source_rows","NCASE_NCON_NTOT_missing_or_nonpositive",
      "NTOT_differs_from_NCASE_plus_NCON_by_more_than_1","NEFF_missing_or_nonpositive","NEFF_exceeds_NTOT"),
      count=c(nrow(d),sum(!is.finite(d$NCASE)|d$NCASE<=0|!is.finite(d$NCON)|d$NCON<=0|!is.finite(d$NTOT)|d$NTOT<=0),
      sum(abs(d$NTOT-d$NCASE-d$NCON)>1,na.rm=TRUE),sum(!is.finite(d$NEFF)|d$NEFF<=0),
      sum(d$NEFF>d$NTOT+1,na.rm=TRUE)))
    fwrite(checks,file.path(out,"source_count_consistency.csv"))
    ok <- with(d,is.finite(FREQ)&FREQ>0&FREQ<1&is.finite(SE)&SE>0&
      is.finite(BETA)&is.finite(NEFF)&NEFF>0&is.finite(NCASE)&NCASE>0&
      is.finite(NCON)&NCON>0&is.finite(NTOT)&NTOT>0&
      toupper(A1)%in%c("A","C","G","T")&toupper(A2)%in%c("A","C","G","T")&toupper(A1)!=toupper(A2))
    z <- copy(d[ok])
    z[, MAF := pmin(FREQ,1-FREQ)]
    z[, approx_neff_SE := 4/(2*FREQ*(1-FREQ)*SE^2)]
    z[, pooled_counts_neff := 4*NCASE*NCON/(NCASE+NCON)]
    z[, SE_over_source_NEFF := approx_neff_SE/NEFF]
    z[, source_NEFF_over_pooled := NEFF/pooled_counts_neff]
    subsets <- list(valid_biallelic=seq_len(nrow(z)),
      common_MAF_ge_005=which(z$MAF>=.05),
      common_MAF_ge_005_absZ_le_2=which(z$MAF>=.05 & abs(z$BETA/z$SE)<=2))
    diagnostics <- list();agreement <- list()
    for(nm in names(subsets)) {
      ix <- subsets[[nm]]
      diagnostics[[nm]] <- rbindlist(lapply(c("approx_neff_SE","SE_over_source_NEFF","source_NEFF_over_pooled"),
        function(v) summarize(z[[v]][ix],v,nm)))
      ratio <- z$SE_over_source_NEFF[ix];ratio <- ratio[is.finite(ratio)]
      agreement[[nm]] <- data.table(subset=nm,n=length(ratio),
        fraction_SE_ratio_within_005=if(length(ratio)) mean(abs(ratio-1)<=.05) else NA_real_,
        fraction_SE_ratio_within_010=if(length(ratio)) mean(abs(ratio-1)<=.10) else NA_real_,
        fraction_SE_ratio_within_020=if(length(ratio)) mean(abs(ratio-1)<=.20) else NA_real_)
    }
    fwrite(rbindlist(diagnostics),file.path(out,"NEFF_formula_diagnostics.csv"))
    fwrite(rbindlist(agreement),file.path(out,"SE_approximation_agreement.csv"))
    chromosome <- z[MAF>=.05 & abs(BETA/SE)<=2,.(n=.N,
      median_SE_over_source_NEFF=median(SE_over_source_NEFF),
      q05_SE_over_source_NEFF=unname(quantile(SE_over_source_NEFF,.05)),
      q95_SE_over_source_NEFF=unname(quantile(SE_over_source_NEFF,.95))),by=CHR]
    setorder(chromosome,CHR)
    fwrite(chromosome,file.path(out,"common_null_SNP_chromosome_diagnostics.csv"))
    targets <- data.table(protein=c("AKT3","CD101","CD40","CGREF1","FES","FURIN","KHK","SIRPA","SNX18","UBE2L6"),
      ID=c("rs4515770","rs12130298","rs4810485","rs2384572","rs1894401","rs2071410","rs2304681","rs6136377","rs2548611","rs28362950"))
    candidate <- merge(targets,z,by="ID",all.x=TRUE,sort=FALSE)
    fwrite(candidate,file.path(out,"candidate_source_diagnostics.csv"))
    patterns <- d[,.(SNP_rows=.N),by=.(NCASE,NCON,NEFF,NTOT)][order(-SNP_rows)]
    fwrite(head(patterns,100L),file.path(out,"top100_source_N_patterns.csv"))
    cat("Distinct per-SNP N patterns: ",nrow(patterns),"\n",sep="")
    print(field_summary)
    print(rbindlist(diagnostics)[quantity=="SE_over_source_NEFF"])

    cat("STEP 3/3: Save the audit interpretation limits\n")
    writeLines(c(
      "SOURCE AUDIT ONLY; no MR/LDSC/SEM/COJO rerun and no input or result overwrite.",
      "Source is Figshare ptsd2024 version 3, verified MD5 bf0fb5b16bf91e640af14784e9d3a829.",
      "SE approximation: 4/[2*f*(1-f)*SE(logOR)^2]. This is diagnostic, not an exact identity.",
      "No capping, truncation to a chosen N, or substitution of source NEFF has been applied.",
      "Common variants (MAF>=0.05) and a common abs(Z)<=2 subset provide additional descriptive checks.",
      "Agreement cannot prove how a provider constructed NEFF; disagreement can reflect frequency/SE rounding, model adjustments, or other meta-analysis details.",
      "Pooled-count effective N is a comparison only. It is not the cohort-summed NEFF required by Genomic SEM.",
      "Public freeze3 code sums cohort-provided neff over cohorts contributing to each SNP, but the exact cohort NEF table and mapping to this public release remain unavailable in the inspected sources.",
      "Supplementary Table 2 reproduces the primary mixed-phenotype counts; it is not verified as the case/control-release cohort accounting table.",
      "Per-SNP N maxima do not establish a unique, de-duplicated participant total.",
      "Do not automatically replace current primary results, source-NEFF sensitivity, or manuscript sample counts based on this diagnostic alone."
    ),file.path(out,"INTERPRETATION_LIMITS.txt"))
    writeLines(c("AUDIT COMPLETE",format(Sys.time())),file.path(out,"SUCCESS.txt"))
    cat("AUDIT COMPLETE. Send the small SEND_BACK folder for review.\n")
  },error=function(e) {
    writeLines(c(format(Sys.time()),conditionMessage(e)),file.path(out,"FAILED.txt"))
    cat("AUDIT STOPPED: ",conditionMessage(e),"\n",sep="")
    stop(e)
  })
  invisible(out)
}

if (!isTRUE(getOption("ptsd.neff_source_audit.load_only", FALSE))) {
  P <- getOption("ptsd.current.paths")
  if (is.null(P)) ptsd_neff_source_audit() else
    ptsd_neff_source_audit(source_file=P$source_casecontrol, output_parent=P$source_input_dir)
}
