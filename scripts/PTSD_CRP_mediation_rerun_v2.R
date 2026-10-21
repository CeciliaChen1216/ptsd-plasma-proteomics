# Current package: path configuration only; scientific calculations preserved.
P <- getOption("ptsd.current.paths")
if (is.null(P)) stop("First source 00_START.R; then use ptsd_run().")
# Exploratory protein -> CRP -> PTSD product-of-coefficients audit, 2026-10-08.
# Does NOT rerun CRP instrument selection or overwrite earlier analyses.
(function() {
  pqtl_path <- P$pqtl_rds
  crp_path <- P$crp_gwas
  step2_path <- P$crp_step2_results
  total_path <- P$main_mr_rds
  total_outcome <- "PTSD_freeze3"
  result_base <- P$result_dir
  pkgs <- c("data.table","R.utils")
  missing <- pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing)) stop("Install packages first: ",paste(missing,collapse=", "))
  library(data.table)
  paths <- c(pqtl_path,crp_path,step2_path,total_path)
  if(any(!file.exists(paths))) stop("Missing input: ",paste(paths[!file.exists(paths)],collapse="; "))
  out <- tempfile(paste0("mediation_rerun_",format(Sys.time(),"%Y%m%d_%H%M%S"),"_"),tmpdir=result_base)
  dir.create(out,recursive=TRUE)
  cat("Output directory: ",out,"\n",sep="")
  con <- file(file.path(out,"RUN.log"),"wt")
  sink(con,split=TRUE);sink(con,type="message")
  on.exit({sink(type="message");sink();close(con)},add=TRUE)
  sf<-tryCatch(sys.frame(1)$ofile,error=function(e)NULL)
  if(!is.null(sf)&&file.exists(sf)) file.copy(sf,file.path(out,"script_snapshot.R.txt"))
  writeLines(c(
    "Exploratory reanalysis of the original ten candidate proteins; no new selection on significance.",
    "First step: single-sentinel Wald ratio with exact rsID and full allele-pair alignment; no automatic strand rescue.",
    "All palindromic A/T or C/G instruments excluded here: supplied CRP data lack effect allele frequency for orientation checking.",
    "Exposure frequency and frequency from another protein assay do not resolve CRP strand ambiguity.",
    "Second step: read beta and multiplicative random-effects SE from the validated cis500kb and cis100kb rerun CSV; no hardcoded estimates.",
    "Indirect estimate=a*b; delta SE=sqrt(b^2*SEa^2+a^2*SEb^2), assuming Cov(a,b)=0.",
    "Zero covariance is an approximation, not verified: CRP data are reused and cross-study sample overlap is unresolved.",
    "No adjustment of CRP->PTSD for each protein using MVMR; products are exploratory indirect-effect estimates, not identified natural causal mediation effects.",
    "PTSD binary outcome is on log-odds scale; causal decomposition and proportion interpretation require additional assumptions.",
    "Signed indirect/total ratio retained only as a descriptive quantity; no absolute value, clipping, percentage-mediated claim or ratio CI.",
    "Negative ratios indicate opposite signs, not negative biological mediation or confirmed suppression.",
    "Normal-approximation 95% CIs and two-sided delta-method P; BH across original 10 candidates separately within each CRP window, n=10.",
    "No IL6 analysis and no network computation in this script.",
    "Method reference: Carter et al. 2021, doi:10.1007/s10654-021-00757-1."
  ),file.path(out,"parameters.txt"))
  tryCatch({
    cat("STEP 1/4: Input and second-step verification\n")
    fi<-file.info(paths)
    fwrite(data.table(path=paths,bytes=fi$size,modified=as.character(fi$mtime)),file.path(out,"input_manifest.csv"))
    print(tools::md5sum(c(pqtl_path,step2_path,total_path)))
    pq<-as.data.table(readRDS(pqtl_path));tot<-as.data.table(readRDS(total_path));b<-fread(step2_path)
    req<-c("protein","SNP","effect_allele","other_allele","beta","se")
    if(!all(req %in% names(pq))) stop("pQTL columns missing: ",paste(setdiff(req,names(pq)),collapse=", "))
    if(!all(c("protein","outcome","mr_beta") %in% names(tot))) stop("Main MR results need protein/outcome/mr_beta columns; no legacy fallback used")
    stopifnot(all(c("analysis","role","beta","se","p","method") %in% names(b)))
    b<-b[analysis %in% c("cis500kb","cis100kb")]
    if(nrow(b)!=2L||anyDuplicated(b$analysis)) stop("Expected exactly one cis500kb and one cis100kb result")
    if(any(!is.finite(b$beta)|!is.finite(b$se)|b$se<=0)) stop("Invalid second-step estimates")
    if(any(abs(2*pnorm(-abs(b$beta/b$se))-b$p)>1e-6)) stop("Second-step P does not match beta/SE")
    fwrite(b,file.path(out,"CRP_to_PTSD_inputs.csv"))
    snps<-c(AKT3="rs4515770",CD40="rs4810485",CGREF1="rs2384572",FES="rs1894401",FURIN="rs2071410",SIRPA="rs6136377",CD101="rs12130298",KHK="rs2304681",SNX18="rs2548611",UBE2L6="rs28362950")
    cat("STEP 2/4: Read CRP and align complete allele pairs\n")
    crp<-fread(crp_path,select=c("variant_id","chromosome","base_pair_location","effect_allele","other_allele","beta","standard_error","p_value"))
    crp<-crp[variant_id %in% unname(snps)]
    fwrite(crp,file.path(out,"CRP_sentinel_raw_rows.csv"))
    audit<-list();first<-list();total_rows<-list()
    for(g in names(snps)) {
      x<-pq[protein==g & SNP==snps[[g]]]
      if(nrow(x)!=1L) stop("Expected exactly one fixed instrument for ",g)
      tt<-tot[protein==g & outcome==total_outcome]
      if(nrow(tt)!=1L) stop("Expected one main PTSD result for ",g," / ",total_outcome)
      if("SNP" %in% names(tt) && as.character(tt$SNP)!=snps[[g]]) stop("Main total effect uses a different SNP for ",g)
      if(!is.finite(tt$mr_beta)) stop("Invalid total effect for ",g)
      total_rows[[g]]<-tt
      info<-data.table(protein=g,SNP=snps[[g]],exposure_EA=toupper(x$effect_allele),exposure_OA=toupper(x$other_allele),bx=x$beta,sx=x$se,total_beta=tt$mr_beta)
      est<-tryCatch({
        if(!is.finite(x$beta)||x$beta==0||!is.finite(x$se)||x$se<=0) stop("Invalid exposure beta/SE")
        y<-crp[variant_id==snps[[g]]]
        if(nrow(y)!=1L) stop("Missing or duplicate CRP SNP ID")
        ea<-info$exposure_EA;oa<-info$exposure_OA;ya<-toupper(y$effect_allele);yo<-toupper(y$other_allele)
        info[,`:=`(CRP_EA=ya,CRP_OA=yo,CRP_beta_raw=y$beta,CRP_SE=y$standard_error,CRP_P_reported=y$p_value)]
        if(!(ea %in% c("A","C","G","T")&&oa %in% c("A","C","G","T")&&ea!=oa)) stop("Non-biallelic instrument")
        if(!((ea==ya&&oa==yo)||(ea==yo&&oa==ya))) stop("Full allele-pair mismatch")
        if(paste(sort(c(ea,oa)),collapse="/") %in% c("A/T","C/G")) stop("Palindromic SNP; CRP EAF unavailable")
        if(!is.finite(y$beta)||!is.finite(y$standard_error)||y$standard_error<=0) stop("Invalid CRP beta/SE")
        sign<-if(ea==ya) 1 else -1
        info[,`:=`(alignment_sign=sign,CRP_beta_aligned=sign*y$beta)]
        aa<-sign*y$beta/x$beta;sa<-y$standard_error/abs(x$beta)
        data.table(protein=g,SNP=snps[[g]],a=aa,se_a=sa,p_a=2*pnorm(-abs(aa/sa)),F_stat=(x$beta/x$se)^2,total_beta=tt$mr_beta,status="testable",reason="")
      },error=function(e)data.table(protein=g,SNP=snps[[g]],a=NA_real_,se_a=NA_real_,p_a=NA_real_,F_stat=(x$beta/x$se)^2,total_beta=tt$mr_beta,status="not_testable",reason=conditionMessage(e)))
      info[,`:=`(status=est$status,reason=est$reason)]
      first[[g]]<-est;audit[[g]]<-info
    }
    a<-rbindlist(first,fill=TRUE)
    fwrite(rbindlist(audit,fill=TRUE),file.path(out,"mediation_harmonization_audit.csv"))
    fwrite(rbindlist(total_rows,fill=TRUE),file.path(out,"main_MR_total_effect_inputs.csv"))
    fwrite(a,file.path(out,"protein_to_CRP_results.csv"))
    cat("STEP 3/4: Exploratory indirect estimates for both CRP windows\n")
    res<-rbindlist(lapply(seq_len(nrow(b)),function(i) {
      d<-copy(a)
      d[,`:=`(analysis=b$analysis[i],role=b$role[i],b=b$beta[i],se_b=b$se[i],p_b=b$p[i])]
      d[,`:=`(indirect_beta=a*b,indirect_se=sqrt(b^2*se_a^2+a^2*se_b^2))]
      d[,`:=`(indirect_lower=indirect_beta-qnorm(.975)*indirect_se,indirect_upper=indirect_beta+qnorm(.975)*indirect_se,indirect_p=2*pnorm(-abs(indirect_beta/indirect_se)))]
      d[,indirect_fdr:=p.adjust(indirect_p,"BH",n=10)]
      d[,signed_indirect_total_ratio:=ifelse(is.finite(total_beta)&total_beta!=0,indirect_beta/total_beta,NA_real_)]
      d[,ratio_note:=ifelse(status!="testable","not_estimable",ifelse(!is.finite(signed_indirect_total_ratio),"total_zero_or_missing",ifelse(signed_indirect_total_ratio<0,"opposite_sign_descriptive_only","descriptive_only_not_mediation_proportion")))]
      d
    }),fill=TRUE)
    fwrite(res,file.path(out,"mediation_CRP_results.csv"))
    cat("STEP 4/4: Save review outputs\n")
    summary<-res[,.(candidate_proteins=.N,testable=sum(status=="testable"),nominal_indirect_p005=sum(indirect_p<.05,na.rm=TRUE),indirect_FDR005=sum(indirect_fdr<.05,na.rm=TRUE)),by=.(analysis,role)]
    fwrite(summary,file.path(out,"summary_counts.csv"));print(summary)
    print(res[role=="primary",.(protein,a,b,indirect_beta,indirect_se,indirect_p,signed_indirect_total_ratio,status,reason)])
    writeLines("Computation complete; review assumptions and harmonization. No confirmatory causal mediation claim.",file.path(out,"COMPLETED.txt"))
    cat("Send mediation_CRP_results.csv, mediation_harmonization_audit.csv, main_MR_total_effect_inputs.csv, summary_counts.csv and RUN.log.\n")
    cat("Output directory: ",out,"\n",sep="")
  },error=function(e){writeLines(conditionMessage(e),file.path(out,"FAILED.txt"));message(conditionMessage(e))},finally={capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))})
})()
