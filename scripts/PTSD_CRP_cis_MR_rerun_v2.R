# Current package: path configuration only; scientific calculations preserved.
P <- getOption("ptsd.current.paths")
if (is.null(P)) stop("First source 00_START.R; then use ptsd_run().")
# CRP -> PTSD cis-MR reanalysis, 2026-10-07.
# This reconstructs instruments under explicitly recorded NEW audit choices;
# it does not claim reproduction of the legacy undocumented instrument list.
# Primary: CRP gene +/-500kb; sensitivity: +/-100kb; P<5e-8, F>10,
# LD r2<0.001 within 10Mb, unphased LD, European reference.
# CRP GRCh37 gene: chr1:159682079-159684379 (Ensembl ENSG00000132693).
# https://www.ensembl.org/feature-explorer/3704ceb1-948d-11ec-a39d-005056b38ce3/gene%3AENSG00000132693
# PLINK2 clump: https://www.cog-genomics.org/plink/2.0/postproc
# No prior results are overwritten; no packages are installed or updated.
# Run: source("path/to/code/scripts/PTSD_CRP_cis_MR_rerun_v2.R", encoding="UTF-8")

run_crp_cis <- function() {
  plink <- P$plink2
  ref <- P$ld_prefix
  crp_file <- P$crp_gwas
  ptsd_file <- P$ptsd_ma
  parent <- P$result_dir
  if (!dir.exists(parent)) stop("Results directory missing")
  out <- tempfile(paste0("crp_cis_rerun_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir=parent)
  if (!dir.create(out)) stop("Cannot create output directory")
  depth <- sink.number()
  con <- file(file.path(out,"RUN.log"), "wt")
  sink(con, split=TRUE); sink(con, type="message")
  on.exit({sink(type="message"); while(sink.number()>depth) sink(); close(con)
           cat("\nOutput directory: ",out,"\n",sep="")}, add=TRUE)
  tryCatch({
    if (!requireNamespace("data.table",quietly=TRUE) || !requireNamespace("R.utils",quietly=TRUE)) stop("data.table and R.utils must be installed")
    suppressPackageStartupMessages(library(data.table))
    writeLines(capture.output(dput(run_crp_cis)),file.path(out,"script_snapshot.R.txt"))
    writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
    paths <- c(plink,crp_file,ptsd_file,paste0(ref,c(".bed",".bim",".fam")))
    if (any(!file.exists(paths))) stop("Missing files: ",paste(paths[!file.exists(paths)],collapse="\n"))
    cat("STEP 1/4: Input verification\n")
    manifest <- data.frame(path=paths,bytes=file.info(paths)$size,md5=unname(tools::md5sum(paths)))
    write.csv(manifest,file.path(out,"input_manifest.csv"),row.names=FALSE)
    ver <- system2(plink,"--version",stdout=TRUE,stderr=TRUE)
    writeLines(ver,file.path(out,"plink_version.txt"));cat(paste(ver,collapse="\n"),"\n")
    params <- list(CRP_accession="GCST90029070",build="GRCh37",gene_chr=1,
       gene_start=159682079L,gene_end=159684379L,primary_flank_bp=500000L,
       sensitivity_flank_bp=100000L,p_threshold=5e-8,F_threshold=10,
       clump_r2=0.001,clump_kb=10000,LD_type="unphased",
       P_handling="Original text parsed to -log10(P) without numeric underflow",
       palindromic="excluded: CRP summary file lacks EAF",
       method="Wald for one SNP; fixed IVW and multiplicative random-effects IVW for >=2",
       outcome_scale="PTSD .ma b assumed log odds, as in main study",
       provenance="New audit reanalysis; original cis-CRP instrument list unavailable")
    dput(params,file=file.path(out,"parameters.txt"))
    bim <- fread(paste0(ref,".bim"),header=FALSE)
    if (ncol(bim)!=6) stop("Unexpected BIM format")
    setnames(bim,c("CHR","SNP","CM","BP","B1","B2"))
    nref <- nrow(fread(paste0(ref,".fam"),header=FALSE))
    cat("LD reference participants:",nref,"\n")
    bim <- bim[CHR==1]
    # Force character input: scientific P values can be smaller than 1e-308.
    crp <- fread(crp_file, colClasses=c(p_value="character"))
    parse_logp <- function(x) {
      x <- tolower(trimws(as.character(x)))
      ans <- rep(NA_real_,length(x))
      sci <- !is.na(x) & grepl("^[+]?[0-9]*[.]?[0-9]+e[+-]?[0-9]+$",x)
      if (any(sci)) {
        mantissa <- suppressWarnings(as.numeric(sub("e.*$","",x[sci])))
        exponent <- suppressWarnings(as.numeric(sub("^.*e","",x[sci])))
        ans[sci] <- -log10(mantissa)-exponent
      }
      ordinary <- !is.na(x) & !sci & grepl("^[+]?[0-9]*[.]?[0-9]+$",x)
      ans[ordinary] <- -log10(suppressWarnings(as.numeric(x[ordinary])))
      ans[!is.finite(ans)|ans<0] <- NA_real_
      ans
    }
    stopifnot(abs(parse_logp("8.69e-1352")-(1352-log10(8.69)))<1e-10,
              abs(parse_logp("5e-8")-(-log10(5e-8)))<1e-12,
              is.na(parse_logp("0")), is.na(parse_logp("<1e-300")))
    required <- c("variant_id","chromosome","base_pair_location","effect_allele","other_allele","beta","standard_error","p_value")
    if (!all(required %in% names(crp))) stop("Unexpected CRP columns")
    e <- crp[chromosome==1 & base_pair_location>=159182079 & base_pair_location<=160184379,
        .(SNP=variant_id,BP=base_pair_location,EA=toupper(effect_allele),OA=toupper(other_allele),bx=beta,sx=standard_error,P_original=p_value)]
    rm(crp);gc()
    e[, SNP:=trimws(as.character(SNP))]
    e[, LOG10P:=parse_logp(P_original)]
    fwrite(e[!is.finite(LOG10P)],file.path(out,"unparsed_P_values.csv"))
    cat("Unparsed/zero/invalid P values:",sum(!is.finite(e$LOG10P)),"\n")
    e[,reason:="eligible"]
    dup <- duplicated(e$SNP)|duplicated(e$SNP,fromLast=TRUE)
    e[dup|is.na(SNP)|!nzchar(SNP),reason:="Missing or duplicate SNP ID"]
    e[reason=="eligible" & (!is.finite(bx)|!is.finite(sx)|sx<=0|!is.finite(LOG10P)|LOG10P<=-log10(5e-8)),reason:="Invalid effect/SE/P or P>=5e-8"]
    e[,F_stat:=(bx/sx)^2]
    e[reason=="eligible" & (!is.finite(F_stat)|F_stat<=10),reason:="F<=10"]
    e[reason=="eligible" & (!EA %in% c("A","C","G","T")|!OA %in% c("A","C","G","T")|EA==OA),reason:="Not a biallelic SNP"]
    e[reason=="eligible" & paste0(EA,OA) %in% c("AT","TA","CG","GC"),reason:="Palindromic; no CRP EAF"]
    # LD filtering requires valid rsID, hg19 position and compatible alleles.
    bdup <- bim$SNP[duplicated(bim$SNP)|duplicated(bim$SNP,fromLast=TRUE)]
    b <- bim[!SNP %in% bdup]
    e <- merge(e,b[,.(SNP,refBP=BP,B1=toupper(B1),B2=toupper(B2))],by="SNP",all.x=TRUE)
    comp <- function(x) chartr("ACGT","TGCA",x)
    e[,ref_ok:= (EA==B1 & OA==B2)|(EA==B2 & OA==B1)|
                  (comp(EA)==B1 & comp(OA)==B2)|(comp(EA)==B2 & comp(OA)==B1)]
    e[reason=="eligible" & (is.na(refBP)|is.na(ref_ok)|!ref_ok|BP!=refBP),reason:="LD reference absent/ambiguous or position/allele mismatch"]
    fwrite(e,file.path(out,"CRP_region_audit.csv"))
    elig <- e[reason=="eligible"]
    if (!nrow(elig)) stop("No eligible cis-CRP instruments before LD clumping")
    cat("Eligible variants before clumping:",nrow(elig),"\n")
    cat("STEP 2/4: Read PTSD effects\n")
    y <- fread(ptsd_file)
    if (!all(c("SNP","A1","A2","b","se") %in% names(y))) stop("Unexpected PTSD .ma columns")
    y <- y[SNP %in% elig$SNP,.(SNP,YA=toupper(A1),YO=toupper(A2),by_raw=b,sy=se)]
    ydup <- y$SNP[duplicated(y$SNP)|duplicated(y$SNP,fromLast=TRUE)]
    y <- y[!SNP %in% ydup]
    # Instruments are clumped using exposure P values before outcome matching.
    cat("STEP 3/4: Clumping and outcome harmonization\n")
    all_results <- list()
    for (flank in c(500000L,100000L)) {
      label <- paste0("cis",flank/1000,"kb")
      d <- elig[BP>=159682079-flank & BP<=159684379+flank]
      if (!nrow(d)) {warning("No eligible variants for ",label);next}
      prefix <- file.path(out,label)
      assoc <- paste0(prefix,"_clump_input.tsv")
      fwrite(d[,.(ID=SNP,LOG10P)],assoc,sep="\t")
      args <- c("--bfile",shQuote(ref),"--chr","1","--clump",shQuote(assoc),
         "--clump-id-field","ID","--clump-p-field","LOG10P",
         "--clump-log10","input-only",
         "--clump-log10-p1",format(-log10(5e-8),digits=16),
         "--clump-log10-p2",format(-log10(5e-8),digits=16),"--clump-r2","0.001","--clump-kb","10000",
         "--clump-unphased","--threads","2","--out",shQuote(prefix))
      writeLines(paste(shQuote(plink),paste(args,collapse=" ")),paste0(prefix,"_command.txt"))
      log <- system2(plink,args,stdout=TRUE,stderr=TRUE)
      writeLines(log,paste0(prefix,"_console.txt"))
      code <- attr(log,"status")
      if (!is.null(code) && code!=0) stop("PLINK failed for ",label,". Read its console/log.")
      cf <- paste0(prefix,".clumps")
      if (!file.exists(cf)) stop("No .clumps file for ",label,". Read PLINK log.")
      cl <- fread(cf)
      if (!"ID" %in% names(cl)) stop("Unexpected .clumps columns")
      selected <- d[SNP %in% cl$ID]
      fwrite(selected,paste0(prefix,"_selected_exposure.csv"))
      h <- merge(selected,y,by="SNP",all.x=TRUE)
      h[,direction:=NA_real_]
      h[(EA==YA & OA==YO)|(comp(EA)==YA & comp(OA)==YO),direction:=1]
      h[(EA==YO & OA==YA)|(comp(EA)==YO & comp(OA)==YA),direction:=-1]
      h[,outcome_status:="included"]
      h[!is.finite(by_raw)|!is.finite(sy)|sy<=0,outcome_status:="Missing/invalid or duplicate PTSD effect"]
      h[outcome_status=="included" & is.na(direction),outcome_status:="PTSD allele mismatch"]
      h[,by_aligned:=by_raw*direction]
      fwrite(h,paste0(prefix,"_harmonization_audit.csv"))
      z <- h[outcome_status=="included"]
      fwrite(z,paste0(prefix,"_instruments_used.csv"))
      if (!nrow(z)) {warning("No usable PTSD instruments for ",label);next}
      z[,`:=`(wald_beta=by_aligned/bx,wald_se=abs(sy/bx))]
      fwrite(z,paste0(prefix,"_SNP_estimates.csv"))
      k <- nrow(z)
      est <- sum(z$bx*z$by_aligned/z$sy^2)/sum(z$bx^2/z$sy^2)
      se_fixed <- sqrt(1/sum(z$bx^2/z$sy^2))
      Q <- if (k>1) sum((z$by_aligned-est*z$bx)^2/z$sy^2) else NA_real_
      se_main <- if(k>1) se_fixed*sqrt(max(1,Q/(k-1))) else se_fixed
      r <- data.table(analysis=label,role=if(flank==500000)"primary" else "sensitivity",
        nsnp=k,selected_before_outcome=nrow(selected),method=if(k==1)"Wald" else "IVW multiplicative random effects",
        beta=est,se=se_main,p=2*pnorm(-abs(est/se_main)),
        OR=exp(est),OR_lower=exp(est-1.96*se_main),OR_upper=exp(est+1.96*se_main),
        se_fixed=se_fixed,p_fixed=2*pnorm(-abs(est/se_fixed)),Q=Q,Q_df=if(k>1)k-1 else NA_integer_,
        Q_p=if(k>1)pchisq(Q,k-1,lower.tail=FALSE) else NA_real_,min_F=min(z$F_stat))
      all_results[[label]] <- r;print(r)
    }
    cat("STEP 4/4: Save results\n")
    if (!length(all_results)) stop("Neither window produced an MR result")
    ans <- rbindlist(all_results)
    fwrite(ans,file.path(out,"CRP_to_PTSD_results.csv"))
    if (!"cis500kb" %in% ans$analysis) stop("Primary window has no estimate; sensitivity must not silently replace it")
    writeLines(c("COMPUTATION COMPLETE; review instruments and warnings before mediation.",
       "Old beta=-0.199 and SE were not used anywhere.",
       "Rounded CRP effects, cis restriction, residual LD, and sample overlap require interpretation."),file.path(out,"SUCCESS.txt"))
    cat("\nCOMPUTATION COMPLETE. Send CRP_to_PTSD_results.csv, cis500kb_instruments_used.csv,\n",
        "cis500kb_harmonization_audit.csv, cis100kb_instruments_used.csv and RUN.log.\n")
  },error=function(e){writeLines(conditionMessage(e),file.path(out,"FAILED.txt"));stop(e)})
  invisible(out)
}
run_crp_cis()
