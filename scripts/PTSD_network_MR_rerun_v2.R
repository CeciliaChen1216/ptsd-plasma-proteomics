# Current package: path configuration only; scientific calculations preserved.
P <- getOption("ptsd.current.paths")
if (is.null(P)) stop("First source 00_START.R; then use ptsd_run().")
# Figure 4B six-pair exploratory network MR audit (2026-10-08).
# Uses existing sentinel instruments; no instrument reselection and no reverse-edge inference.
# All paths are configurable below. No original inputs/results are overwritten.
(function() {
  pqtl_path <- P$pqtl_rds
  raw_base <- P$protein_raw_base
  result_base <- P$result_dir
  legacy_path <- file.path(result_base,"network_mr_results.csv")
  packages <- c("data.table","R.utils")
  missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing)) stop("Install packages first: ",paste(missing,collapse=", "))
  library(data.table)
  if(!file.exists(pqtl_path)) stop("Missing ",pqtl_path)
  dir.create(result_base,recursive=TRUE,showWarnings=FALSE)
  out <- tempfile(paste0("network_rerun_",format(Sys.time(),"%Y%m%d_%H%M%S"),"_"),tmpdir=result_base)
  dir.create(out)
  cat("Output directory: ",out,"\n",sep="")
  con <- file(file.path(out,"RUN.log"),"wt")
  sink(con,split=TRUE);sink(con,type="message")
  on.exit({sink(type="message");sink();close(con)},add=TRUE)
  sf <- tryCatch(sys.frame(1)$ofile,error=function(e)NULL)
  if(!is.null(sf)&&file.exists(sf)) file.copy(sf,file.path(out,"script_snapshot.R.txt"))
  writeLines(c(
    "Six directional pairs from existing Figure 4B; fixed sentinel instruments, no selection on new results.",
    "Discovery files selected by exact protein assay and exposure chromosome; never take the first matching file.",
    "Match GRCh38 GENPOS and complete allele pair; cross-check GRCh37 position embedded in raw ID across exposure/outcome.",
    "BETA refers to ALLELE1. Outcomes aligned to exposure effect_allele. No automatic strand-complement rescue.",
    "Require ADD, INFO>=0.8, positive finite SE/N and biallelic SNVs.",
    "Frequency difference after alignment <=0.10; palindromes additionally require the alternative orientation to be worse by >0.10.",
    "Exposure raw beta/SE must reproduce sentinel RDS within max(1e-8,abs(value)*1e-5).",
    "Wald beta=outcome_beta/exposure_beta; first-order SE=outcome_SE/abs(exposure_beta).",
    "BH and Bonferroni correction across six specified pairs, including untestable pairs in the family size.",
    "Both traits use UKB-PPP discovery data: sample overlap/covariance is not resolved. SE treats exposure association as fixed.",
    "These exploratory estimates alone do not establish directed protein regulation, mediation, or independent causal network edges.",
    "A non-testable pair is reported as such, never substituted by another SNP or protein assay."
  ),file.path(out,"parameters.txt"))
  tryCatch({
    cat("STEP 1/4: Fixed pairs, instruments and file verification\n")
    pairs <- data.table(upstream=c("AKT3","CD40","FES","CD40","SIRPA","FURIN"),
                        downstream=c("FURIN","UBE2L6","UBE2L6","FES","CD40","CD40"))
    sentinel <- c(AKT3="rs4515770",CD40="rs4810485",FES="rs1894401",SIRPA="rs6136377",FURIN="rs2071410")
    assays <- c(AKT3="AKT3_Q9Y243_OID21197_v1_Oncology",CD40="CD40_P25942_OID20724_v1_Inflammation",
                FES="FES_P07332_OID21207_v1_Oncology",FURIN="FURIN_P09958_OID21514_v1_Oncology",
                SIRPA="SIRPA_P78324_OID20304_v1_Cardiometabolic",UBE2L6="UBE2L6_O14933_OID30321_v1_Cardiometabolic_II")
    pqtl <- as.data.table(readRDS(pqtl_path))
    req <- c("protein","SNP","chr","pos_hg38","effect_allele","other_allele","beta","se","eaf")
    if(!all(req %in% names(pqtl))) stop("Missing sentinel columns: ",paste(setdiff(req,names(pqtl)),collapse=", "))
    inst <- rbindlist(lapply(names(sentinel),function(g) {
      d <- pqtl[protein==g & SNP==sentinel[[g]]]
      if(nrow(d)!=1L) stop("Expected exactly one fixed sentinel for ",g)
      d[,..req]
    }))
    inst[,`:=`(effect_allele=toupper(effect_allele),other_allele=toupper(other_allele))]
    if(any(!is.finite(inst$beta)|inst$beta==0|!is.finite(inst$se)|inst$se<=0|!is.finite(inst$eaf)|inst$eaf<=0|inst$eaf>=1)) stop("Invalid sentinel statistics")
    fwrite(inst,file.path(out,"fixed_instruments.csv"))
    requests <- rbindlist(lapply(seq_len(nrow(pairs)),function(i) {
      x <- inst[protein==pairs$upstream[i]]
      data.table(protein=c(pairs$upstream[i],pairs$downstream[i]),chr=x$chr,pos=x$pos_hg38)
    }))
    requests <- unique(requests)
    requests[,path:=file.path(raw_base,protein,paste0("discovery_chr",chr,"_",assays[protein],".gz"))]
    fwrite(requests,file.path(out,"file_requests.csv"))
    if(any(!file.exists(requests$path))) stop("Missing exact files: ",paste(unique(requests$path[!file.exists(requests$path)]),collapse="; "))
    paths <- unique(c(pqtl_path,requests$path)); fi<-file.info(paths)
    manifest <- data.table(path=paths,bytes=fi$size,modified=as.character(fi$mtime),md5=unname(tools::md5sum(paths)))
    fwrite(manifest,file.path(out,"input_manifest.csv"))
    cat("STEP 2/4: Extract all requested positions from exact chromosome files\n")
    rows <- list()
    for(f in unique(requests$path)) {
      pos_needed <- requests[path==f,pos]
      d <- fread(f,select=c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","TEST","BETA","SE"))
      d <- d[GENPOS %in% pos_needed]
      d[,`:=`(source_file=f,ALLELE0=toupper(ALLELE0),ALLELE1=toupper(ALLELE1))]
      rows[[f]] <- d
    }
    raw <- rbindlist(rows,use.names=TRUE)
    fwrite(raw,file.path(out,"extracted_raw_rows.csv"))
    find_row <- function(g,x) {
      f <- requests[protein==g & chr==x$chr & pos==x$pos_hg38,path]
      if(length(f)!=1L) stop("Ambiguous file request for ",g)
      q <- raw[source_file==f & CHROM==x$chr & GENPOS==x$pos_hg38 & TEST=="ADD"]
      ea<-x$effect_allele;oa<-x$other_allele
      if(!(ea %in% c("A","C","G","T") && oa %in% c("A","C","G","T") && ea!=oa)) stop("Non-biallelic sentinel")
      q <- q[(ALLELE1==ea & ALLELE0==oa)|(ALLELE1==oa & ALLELE0==ea)]
      if(nrow(q)!=1L) stop("Expected one full allele-pair match in ",g,"; found ",nrow(q))
      if(!is.finite(q$BETA)||!is.finite(q$SE)||q$SE<=0||!is.finite(q$INFO)||q$INFO<.8||!is.finite(q$N)||q$N<=0||!is.finite(q$A1FREQ)||q$A1FREQ<=0||q$A1FREQ>=1) stop("Invalid raw statistics/INFO in ",g)
      sign <- if(q$ALLELE1==ea) 1 else -1
      q[,`:=`(aligned_beta=sign*BETA,aligned_eaf=if(sign==1) A1FREQ else 1-A1FREQ,alignment_sign=sign)]
      q
    }
    cat("STEP 3/4: Harmonization and six Wald estimates\n")
    ans <- list();audit<-list()
    for(i in seq_len(nrow(pairs))) {
      a<-pairs$upstream[i];b<-pairs$downstream[i];x<-inst[protein==a]
      info<-data.table(upstream=a,downstream=b,SNP=x$SNP,exposure_EA=x$effect_allele,exposure_OA=x$other_allele)
      result <- tryCatch({
        ex<-find_row(a,x);oy<-find_row(b,x)
        id1<-strsplit(ex$ID,":",fixed=TRUE)[[1]];id2<-strsplit(oy$ID,":",fixed=TRUE)[[1]]
        if(length(id1)<4||length(id2)<4||!identical(id1[1:2],id2[1:2])) stop("Embedded chromosome/GRCh37 ID mismatch")
        info[,`:=`(pos_hg38=x$pos_hg38,raw_exposure_ID=ex$ID,raw_outcome_ID=oy$ID,
                   bx=x$beta,sx=x$se,raw_bx=ex$aligned_beta,raw_sx=ex$SE,
                   by_raw=oy$BETA,by_aligned=oy$aligned_beta,sy=oy$SE,
                   exposure_eaf=x$eaf,raw_exposure_eaf=ex$aligned_eaf,outcome_eaf_aligned=oy$aligned_eaf,
                   outcome_alignment_sign=oy$alignment_sign,N_exposure=ex$N,N_outcome=oy$N)]
        if(abs(ex$aligned_beta-x$beta)>max(1e-8,abs(x$beta)*1e-5)||abs(ex$SE-x$se)>max(1e-8,abs(x$se)*1e-5)) stop("Raw exposure beta/SE do not reproduce fixed instrument")
        if(abs(ex$aligned_eaf-x$eaf)>.01) stop("Raw exposure frequency differs from instrument")
        direct<-abs(oy$aligned_eaf-ex$aligned_eaf);opposite<-abs(1-oy$aligned_eaf-ex$aligned_eaf)
        pal<-paste(sort(c(x$effect_allele,x$other_allele)),collapse="/") %in% c("A/T","C/G")
        info[,`:=`(palindromic=pal,frequency_distance=direct,opposite_frequency_distance=opposite)]
        if(direct>.10) stop("Aligned allele frequency difference >0.10")
        if(pal && opposite-direct<=.10) stop("Palindromic orientation ambiguous by frequency")
        beta<-oy$aligned_beta/x$beta;se<-oy$SE/abs(x$beta);z<-beta/se
        data.table(upstream=a,downstream=b,SNP=x$SNP,beta=beta,se=se,z=z,
                   pval=2*pnorm(-abs(z)),log10p=-(log(2)+pnorm(-abs(z),log.p=TRUE))/log(10),
                   lower=beta-qnorm(.975)*se,upper=beta+qnorm(.975)*se,F_stat=(x$beta/x$se)^2,
                   status="testable",reason="")
      },error=function(e)data.table(upstream=a,downstream=b,SNP=x$SNP,beta=NA_real_,se=NA_real_,z=NA_real_,
               pval=NA_real_,log10p=NA_real_,lower=NA_real_,upper=NA_real_,F_stat=(x$beta/x$se)^2,
               status="not_testable",reason=conditionMessage(e)))
      info[,`:=`(status=result$status,reason=result$reason)]
      audit[[i]]<-info;ans[[i]]<-result
      cat(a," -> ",b,": ",result$status," ",result$reason,"\n",sep="")
    }
    net<-rbindlist(ans,fill=TRUE)
    net[,`:=`(fdr=p.adjust(pval,"BH",n=6),p_bonferroni=p.adjust(pval,"bonferroni",n=6))]
    fwrite(net,file.path(out,"network_mr_results_recomputed.csv"))
    fwrite(rbindlist(audit,fill=TRUE),file.path(out,"network_harmonization_audit.csv"))
    cat("STEP 4/4: Legacy comparison and output\n")
    if(file.exists(legacy_path)) {
      old<-fread(legacy_path)
      if(all(c("upstream","downstream","beta","se","pval") %in% names(old)) && !anyDuplicated(old[,.(upstream,downstream)])) {
        old<-old[,.(upstream,downstream,beta_legacy=beta,se_legacy=se,pval_legacy=pval)]
        cmp<-merge(net,old,by=c("upstream","downstream"),all.x=TRUE,sort=FALSE)
        cmp[,`:=`(beta_difference=beta-beta_legacy,se_difference=se-se_legacy,pval_difference=pval-pval_legacy)]
        fwrite(cmp,file.path(out,"network_mr_results_comparison.csv"))
      } else writeLines("Legacy schema missing required columns or duplicate edges; compare manually.",file.path(out,"legacy_comparison_note.txt"))
    } else writeLines(paste("Legacy result absent:",legacy_path),file.path(out,"legacy_comparison_note.txt"))
    print(net)
    writeLines(paste("Audit completed; testable",sum(net$status=="testable"),"of 6. Review all non-testable pairs and overlap limitation."),file.path(out,"COMPLETED.txt"))
    cat("Send network_mr_results_recomputed.csv, network_harmonization_audit.csv, network_mr_results_comparison.csv (if present), and RUN.log.\n")
    cat("Output directory: ",out,"\n",sep="")
  },error=function(e){writeLines(conditionMessage(e),file.path(out,"FAILED.txt"));message(conditionMessage(e))},
  finally={capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))})
})()
