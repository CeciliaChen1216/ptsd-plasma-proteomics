# PTSD follow-up, 2026-10-16. One entry file; original inputs/results are read-only.
# Maintenance fix 2026-10-17: participant overlap uses sample_id, not the reserved
# data.table constructor argument key. No SEM calculations were changed.
# Maintenance fix 2026-10-17: explicitly intersect PLINK2 coordinate and ID filters;
# include extraction/frequency logs and complete SEM result tables in SEND_BACK.
# Run in a restarted R session:
# Use the bundle: source("path/to/code/00_START.R", encoding="UTF-8")
# In this bundle use 00_START.R and ptsd_run("sem_neff") / ptsd_run("cojo").
# No package installation, downloads, manuscript changes, or repository uploads.
# Outputs go to a NEW configured results directory/sem_cojo_followup_<timestamp>_<random>/.
#
# Scope: PTSD source-NEFF SEM sensitivity + PTSD-only COJO reference sensitivity.
# NEFF is NOT assumed proven to be sum of cohort-specific effective sample sizes.
# MDD/ANX N, prevalence and model choices are retained for a controlled comparison.
# COJO with 503 reference participants is exploratory; GCTA recommends larger
# reference samples. No conditional-Wald protein causal estimates are calculated.
# Trauma COJO remains pending: its source/frequency handling is not fixed here.
#
# Optional, only if a module needs rerunning later:
# options(ptsd.followup.mode="cojo")   # default "both"; also "sem"
# options(ptsd.followup.gcta="/path/to/gcta64") # if auto-discovery is ambiguous
# options(ptsd.followup.load_only=TRUE) # define functions without running (testing)
#
# References (retrieved 2026-10-16):
# https://github.com/GenomicSEM/GenomicSEM/wiki/2.1-Calculating-Sum-of-Effective-Sample-Size-and-Preparing-GWAS-Summary-Statistics
# https://github.com/GenomicSEM/GenomicSEM/blob/master/PATCHNOTES.md
# https://yanglab.westlake.edu.cn/software/gcta/#COJO
# Official source: https://doi.org/10.6084/m9.figshare.26349322.v3

ptsd_followup_note <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), " | ", ..., "\n", sep="")
  flush.console()
}

ptsd_read_header <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con))
  lines <- character()
  repeat {
    z <- readLines(con, n=1L, warn=FALSE)
    if (!length(z)) stop("No #CHROM header in PTSD source")
    lines <- c(lines,z)
    if (startsWith(z,"#CHROM")) break
    if (length(lines)>10000L) stop("Unexpected PTSD header")
  }
  lines
}

ptsd_summary_n <- function(x, field) {
  v <- as.numeric(x); v <- v[is.finite(v)]
  if (!length(v)) stop("No valid values for ",field)
  data.table::data.table(field=field, n_valid=length(v), min=min(v),
      q01=unname(quantile(v,.01)), median=median(v), q99=unname(quantile(v,.99)), max=max(v))
}

ptsd_prepare_followup <- function(P) {
  od <- file.path(P$followup_out,"INPUTS"); dir.create(od)
  raw_path <- P$source_casecontrol
  if (!file.exists(raw_path)) stop("Missing downloaded case-control source: ",raw_path)
  ptsd_followup_note("Checking the downloaded case-control file...")
  md <- unname(tools::md5sum(raw_path))
  write.csv(data.frame(path=raw_path,bytes=file.info(raw_path)$size,md5=md),
            file.path(od,"PTSD_source_manifest.csv"),row.names=FALSE)
  if (!identical(md,"bf0fb5b16bf91e640af14784e9d3a829"))
    stop("Case-control source MD5 differs from the audited official release. No analysis started.")
  header <- ptsd_read_header(raw_path)
  writeLines(header,file.path(od,"PTSD_source_header.txt"))
  if (!any(grepl("pgc-ptsd-2024-european-casecontrol",header,fixed=TRUE)))
    stop("Not the verified case-control release")
  ptsd_followup_note("Reading PTSD source once; this may take several minutes...")
  x <- data.table::fread(raw_path,skip="#CHROM",na.strings=c("NA","."))
  data.table::setnames(x,"#CHROM","CHR")
  need <- c("CHR","ID","POS","A1","A2","FREQ","BETA","SE","P","NCASE","NCON","NEFF","NTOT")
  if (!all(need %in% names(x))) stop("PTSD source has unexpected columns")
  if (nrow(x)!=6545577L) stop("PTSD row count differs from the audited release")
  if (anyNA(x$ID) || any(!nzchar(x$ID)) || anyDuplicated(x$ID)) stop("Missing/duplicate PTSD SNP identifiers")
  data.table::set(x,j="A1",value=toupper(x$A1)); data.table::set(x,j="A2",value=toupper(x$A2))
  if (any(!is.finite(x$NEFF)|x$NEFF<=0|!is.finite(x$NTOT)|x$NTOT<=0))
    stop("Invalid source sample sizes; inspect before using NEFF")
  if (any(!is.finite(x$NCASE)|x$NCASE<=0|!is.finite(x$NCON)|x$NCON<=0))
    stop("Invalid case/control sample sizes")
  if (any(abs(x$NCASE+x$NCON-x$NTOT)>1)) stop("Source NTOT differs from NCASE+NCON")
  if (any(x$NEFF>x$NTOT+1)) stop("Source NEFF exceeds NTOT")
  ns <- data.table::rbindlist(lapply(c("NCASE","NCON","NEFF","NTOT"),function(z) ptsd_summary_n(x[[z]],z)))
  data.table::fwrite(ns,file.path(od,"PTSD_sample_size_summary.csv"))
  pooled <- 4*x$NCASE*x$NCON/x$NTOT
  diag <- data.table::rbindlist(list(ptsd_summary_n(pooled,"pooled_count_NEFF_diagnostic_only"),
             ptsd_summary_n(x$NEFF/pooled,"source_NEFF_divided_by_pooled_NEFF")))
  data.table::fwrite(diag,file.path(od,"NEFF_diagnostics.csv"))
  # Pooled-count NEFF is only a diagnostic and is never used as analysis N.
  writeLines(c(
    "SENSITIVITY ANALYSIS, NOT AN AUTOMATIC CORRECTION TO THE MANUSCRIPT",
    "PTSD SEM: N=source NEFF, sample.prev=0.5; log-OR BETA and SE unchanged.",
    "Separate source diagnostics support the NEFF information scale and a cohort-additive structure; complete absolute cohort weights and unique global N remain unresolved (S24c).",
    "Source header counts and per-SNP counts disagree; neither per-SNP maxima nor their sum establish the unique participant total.",
    "No pooled NCASE/NCON formula, capping, or invented global sample count is substituted.",
    "MDD/ANX sample-size and population-prevalence assumptions are held fixed, not independently revalidated.",
    "COJO uses source NTOT, as in the prior input. It never reads the SEM NEFF input.",
    "The 503-person COJO reference remains small; output is exploratory.",
    "Trauma COJO is not rerun; raw outcome frequencies and provenance need separate review."
  ),file.path(P$followup_out,"INTERPRETATION_LIMITS.txt"))
  anchors <- data.table::data.table(ID=c("rs1894401","rs2071410"),
      expected_POS=c(91429042L,91420940L),expected_A1=c("A","C"),expected_A2=c("G","G"))
  aa <- merge(anchors,x,by="ID",all.x=TRUE)
  data.table::fwrite(aa,file.path(od,"FES_FURIN_source_anchors.csv"))
  if (anyNA(aa$POS) || any(aa$CHR!=15 | aa$POS!=aa$expected_POS |
        aa$A1!=aa$expected_A1 | aa$A2!=aa$expected_A2)) stop("GRCh37 anchor verification failed")
  targets <- c(anchors$ID,"rs4515770","rs4810485","rs2384572","rs2304681",
               "rs6136377","rs12130298","rs2548611","rs28362950","rs641325","rs13387644")
  data.table::fwrite(x[ID %in% targets],file.path(od,"candidate_source_rows.csv"))
  # This controlled comparison must not silently change beta/SE/alleles as well as N.
  if (P$mode %in% c("sem","both")) {
    ptsd_followup_note("Checking all rows against the previously used PTSD .ma file...")
    if (!file.exists(P$ptsd_ma)) stop("Original PTSD .ma file missing; cannot isolate the sample-size change")
    previous <- data.table::fread(P$ptsd_ma)
    fields <- c(SNP="ID",A1="A1",A2="A2",freq="FREQ",b="BETA",se="SE",p="P",N="NTOT")
    if (!all(names(fields)%in%names(previous)) || anyDuplicated(previous$SNP))
      stop("Unexpected original .ma fields or duplicate SNPs")
    ix <- match(previous$SNP,x$ID)
    if (nrow(previous)!=nrow(x) || anyNA(ix)) stop("Original .ma and downloaded source have different SNP sets")
    cmp <- data.table::data.table(field=character(),mismatched_rows=integer())
    for (field in names(fields)) {
      lhs <- previous[[field]]; rhs <- x[[fields[[field]]]][ix]
      ok <- if (field %in% c("SNP","A1","A2")) as.character(lhs)==as.character(rhs)
            else is.finite(lhs)&is.finite(rhs)&
                 abs(lhs-rhs)<=1e-8*pmax(abs(lhs),abs(rhs),.Machine$double.xmin)
      bad <- which(is.na(ok)|!ok)
      cmp <- rbind(cmp,data.table::data.table(field=field,mismatched_rows=length(bad)))
      if(length(bad)) {
        keep <- head(bad,20L)
        data.table::fwrite(data.table::data.table(SNP=previous$SNP[keep],old=lhs[keep],source=rhs[keep]),
               file.path(od,paste0("source_mismatch_",field,".csv")))
      }
    }
    data.table::fwrite(cmp,file.path(od,"PTSD_full_source_comparison.csv"))
    if (any(cmp$mismatched_rows>0L)) stop("Old .ma differs from source; see PTSD_full_source_comparison.csv before interpreting a NEFF-only rerun")
    rm(previous,ix,lhs,rhs,ok); gc()
  }
  lo <- min(anchors$expected_POS)-1000000L; hi <- max(anchors$expected_POS)+1000000L
  region <- x[CHR==15 & POS>=lo & POS<=hi]
  data.table::fwrite(region,file.path(od,"PTSD_chr15_GRCh37.csv"))
  ma_path <- file.path(od,"PTSD_NEFF_SEM_only.ma")
  if (P$mode %in% c("sem","both")) {
    data.table::fwrite(x[,.(SNP=ID,A1,A2,freq=FREQ,b=BETA,se=SE,p=P,N=NEFF)],ma_path,sep="\t")
  }
  rm(x,pooled); gc()
  list(sem_ma=ma_path,region=region,anchors=anchors,lower=lo,upper=hi)
}

ptsd_read_fam <- function(prefix) {
  f <- paste0(prefix,".fam")
  if (!file.exists(f)) return(NULL)
  x <- data.table::fread(f,header=FALSE,colClasses="character")
  if (ncol(x)!=6L) stop("Unexpected FAM layout: ",f)
  data.table::setnames(x,c("FID","IID","PAT","MAT","SEX","PHENO"))
  if (anyDuplicated(x[,paste(FID,IID,sep="\t")])) stop("Duplicate FID/IID in ",f)
  x
}

ptsd_find_gcta <- function(P,od) {
  explicit <- getOption("ptsd.followup.gcta","")
  if (length(explicit)!=1L || is.na(explicit)) stop("Invalid ptsd.followup.gcta option")
  candidates <- character()
  if (nzchar(explicit)) {
    if (!file.exists(explicit)) stop("Explicit GCTA path does not exist: ",explicit)
    candidates <- explicit
  } else {
    candidates <- unname(Sys.which(c("gcta64","gcta64.exe","gcta","gcta.exe")))
    home_dir <- Sys.getenv("USERPROFILE")
    dirs <- unique(c(file.path(P$project_root,"tools"),file.path(P$workspace_root,"Longcovid/tools"),file.path(P$workspace_root,"tools"), file.path(P$workspace_root,"synapse/tools"),
        file.path(P$workspace_root,"reference/tools"),file.path(P$workspace_root,"tools"),if(nzchar(home_dir)) file.path(home_dir,"Downloads")))
    dirs <- dirs[dir.exists(dirs)]
    candidates <- c(candidates,unlist(lapply(dirs,function(d)
      list.files(d,pattern="^gcta.*[.]exe$",recursive=TRUE,full.names=TRUE,ignore.case=TRUE))))
    # Also inspect top-level executables, without scanning GWAS directories recursively.
    roots <- c(P$project_root,file.path(P$workspace_root,"Longcovid"),file.path(P$workspace_root,"synapse"),file.path(P$workspace_root,"reference"))
    candidates <- c(candidates,unlist(lapply(roots[dir.exists(roots)],function(d)
      list.files(d,pattern="^gcta.*[.]exe$",full.names=TRUE,ignore.case=TRUE))))
  }
  candidates <- unique(candidates[nzchar(candidates) & file.exists(candidates)])
  tab <- data.table::data.table(path=candidates,md5=unname(tools::md5sum(candidates)))
  data.table::fwrite(tab,file.path(od,"GCTA_candidates.csv"))
  if (!length(candidates)) stop("GCTA not found in known tool folders. SEM can still run. Send GCTA_candidates.csv and your gcta executable path.")
  distinct <- tab[!duplicated(md5)]
  if (nrow(distinct)>1L) stop("Several different GCTA binaries were found. SEM can still run. Send GCTA_candidates.csv so the intended binary can be selected.")
  normalizePath(distinct$path[1],winslash="/",mustWork=TRUE)
}

ptsd_run_external <- function(exe,args,log_path) {
  cmd <- paste(shQuote(exe),paste(vapply(args,shQuote,character(1)),collapse=" "))
  writeLines(cmd,paste0(log_path,".command.txt"))
  ans <- suppressWarnings(system2(exe,args=vapply(args,shQuote,character(1)),stdout=TRUE,stderr=TRUE))
  writeLines(ans,log_path)
  status <- attr(ans,"status"); if(is.null(status)) status <- 0L
  if (status!=0L) stop("External program failed (",status,"): ",log_path)
  invisible(ans)
}

ptsd_cojo_effect <- function(path,snp,A1,A2) {
  if (!file.exists(path)) return(data.table::data.table(status="conditional_output_missing"))
  z <- data.table::fread(path)
  need <- c("SNP","refA","bC","bC_se","pC")
  if (!all(need %in% names(z))) return(data.table::data.table(status="unexpected_COJO_columns"))
  z <- z[SNP==snp]
  if (nrow(z)!=1L) return(data.table::data.table(status="target_not_estimable_or_excluded"))
  s <- if (toupper(z$refA)==A1) 1 else if (toupper(z$refA)==A2) -1 else NA_real_
  if (!is.finite(s) || !is.finite(z$bC) || !is.finite(z$bC_se) || z$bC_se<=0 ||
      !is.finite(z$pC) || z$pC<0 || z$pC>1)
    return(data.table::data.table(status="invalid_or_unaligned_conditional_effect"))
  data.table::data.table(status="estimated",beta_conditional=s*z$bC,se_conditional=z$bC_se,
                        p_conditional=z$pC,COJO_refA=z$refA)
}

run_ptsd_cojo_followup <- function(P,inp) {
  od <- file.path(P$followup_out,"COJO_PTSD"); dir.create(od)
  current <- ptsd_read_fam(P$ld_prefix)
  old_prefix <- file.path(P$project_root,"results/cojo_chr15/1000G_EUR_chr15_region")
  old <- ptsd_read_fam(old_prefix)
  if (is.null(current)) stop("Current LD FAM missing")
  ns <- data.table::data.table(panel="current",prefix=P$ld_prefix,n=nrow(current),
       founders=sum(current$PAT=="0" & current$MAT=="0"),
       fam_md5=unname(tools::md5sum(paste0(P$ld_prefix,".fam"))))
  if (!is.null(old)) ns <- rbind(ns,data.table::data.table(panel="old_region",prefix=old_prefix,n=nrow(old),
      founders=sum(old$PAT=="0" & old$MAT=="0"),fam_md5=unname(tools::md5sum(paste0(old_prefix,".fam")))))
  data.table::fwrite(ns,file.path(od,"LD_reference_audit.csv"))
  if (!is.null(old)) {
    keys1 <- paste(current$FID,current$IID,sep="/")
    keys2 <- paste(old$FID,old$IID,sep="/")
    sample_ids <- union(keys1,keys2)
    data.table::fwrite(data.table::data.table(sample_id=sample_ids,
      in_current=sample_ids%in%keys1,in_old=sample_ids%in%keys2),
      file.path(od,"LD_participant_overlap.csv"))
  }
  writeLines(c("EXPLORATORY SMALL-REFERENCE COJO SENSITIVITY.",
     "Founder status does not establish genetic unrelatedness.",
     "This output tests locus-level conditional association, not independent protein causality.",
     "Only PTSD is rerun. Prior trauma conditional results are still unverified."),file.path(od,"LIMITATIONS.txt"))
  if (nrow(current)!=503L || ns$fam_md5[1]!="7dc68f97189d131278b1e3c39951dbf4")
    stop("Current reference is not the audited 503-participant FAM")
  required <- c(P$plink2,paste0(P$ld_prefix,c(".bed",".bim",".fam")))
  if(any(!file.exists(required))) stop("Missing LD/PLINK input: ",paste(required[!file.exists(required)],collapse=", "))
  gcta <- ptsd_find_gcta(P,od)
  ptsd_followup_note("COJO: checking GRCh37 positions and allele pairs...")
  b <- data.table::fread(paste0(P$ld_prefix,".bim"),header=FALSE)
  if(ncol(b)!=6L) stop("Unexpected BIM layout")
  data.table::setnames(b,c("CHR_ref","ID","CM","BP_ref","B1","B2"))
  # Remove ALL copies of any duplicate ID, including duplicates outside the region.
  duplicate_ids <- b$ID[duplicated(b$ID)|duplicated(b$ID,fromLast=TRUE)]
  br <- b[as.character(CHR_ref)=="15" & BP_ref>=inp$lower & BP_ref<=inp$upper]
  ba <- b[ID %in% inp$anchors$ID]
  data.table::fwrite(ba,file.path(od,"reference_anchor_rows.csv"))
  if(nrow(ba)!=2L || any(ba$ID %in% duplicate_ids)) stop("Missing or duplicate anchor SNPs in reference")
  ba <- merge(ba,inp$anchors,by="ID")
  if(any(as.character(ba$CHR_ref)!="15" | ba$BP_ref!=ba$expected_POS))
    stop("Reference anchor positions are not verified GRCh37; no build conversion guessed")
  rm(b); gc()
  br <- br[!ID %in% duplicate_ids & !is.na(ID) & nzchar(ID)]
  z <- merge(inp$region,br,by="ID",all.x=TRUE)
  z[, `:=`(B1=toupper(B1),B2=toupper(B2))]
  z[, reason := "included"]
  z[is.na(BP_ref),reason := "absent_or_duplicate_reference_ID"]
  z[reason=="included" & POS!=BP_ref,reason := "GRCh37_position_mismatch"]
  z[reason=="included" & !(A1 %in% c("A","C","G","T") & A2 %in% c("A","C","G","T") & A1!=A2),
    reason := "not_biallelic_SNV"]
  z[reason=="included" & !((A1==B1 & A2==B2)|(A1==B2 & A2==B1)),reason := "allele_pair_mismatch"]
  z[reason=="included" & (!is.finite(BETA)|!is.finite(SE)|SE<=0|!is.finite(P)|P<0|P>1|
      !is.finite(FREQ)|FREQ<=0|FREQ>=1|!is.finite(NTOT)|NTOT<=0),reason := "invalid_GWAS_values"]
  z[, palindromic := paste0(A1,A2)%in%c("AT","TA","CG","GC")]
  z[reason=="included" & palindromic & pmin(FREQ,1-FREQ)>.42,
    reason := "palindromic_near_half_frequency"]
  data.table::fwrite(z,file.path(od,"harmonization_before_frequency.csv"))
  if (!all(inp$anchors$ID %in% z[reason=="included",ID])) stop("Anchor did not pass reference alignment")
  extract <- file.path(od,"eligible_ids.txt"); writeLines(z[reason=="included",ID],extract)
  refout <- file.path(od,"EUR503_chr15")
  ptsd_run_external(P$plink2,c("--bfile",P$ld_prefix,"--chr","15","--from-bp",as.character(inp$lower),
       "--to-bp",as.character(inp$upper),"--extract",extract,"--force-intersect",
       "--make-bed","--out",refout,"--threads","2"),
       file.path(od,"PLINK_extract_console.txt"))
  frout <- file.path(od,"EUR503_frequency")
  ptsd_run_external(P$plink2,c("--bfile",refout,"--freq","--out",frout,"--threads","2"),
                    file.path(od,"PLINK_frequency_console.txt"))
  fq <- data.table::fread(paste0(frout,".afreq"))
  if(!all(c("ID","REF","ALT","ALT_FREQS")%in%names(fq))) stop("Unexpected PLINK2 frequency format")
  if(anyDuplicated(fq$ID)) stop("Duplicate IDs in PLINK frequency output")
  fq <- fq[,.(ID,REF=toupper(REF),ALT=toupper(ALT),ALT_FREQS=as.numeric(ALT_FREQS))]
  z <- merge(z,fq,by="ID",all.x=TRUE)
  z[, ref_EAF := ifelse(A1==ALT & A2==REF,ALT_FREQS,
                        ifelse(A1==REF & A2==ALT,1-ALT_FREQS,NA_real_))]
  z[, frequency_difference := abs(FREQ-ref_EAF)]
  z[reason=="included" & (!is.finite(ref_EAF)|ref_EAF<=0|ref_EAF>=1),reason := "invalid_reference_frequency"]
  z[reason=="included" & frequency_difference>.10,reason := "EAF_difference_above_0.10"]
  data.table::fwrite(z,file.path(od,"harmonization_audit.csv"))
  data.table::fwrite(z[ID %in% inp$anchors$ID],file.path(od,"anchor_frequency_check.csv"))
  if(!all(inp$anchors$ID %in% z[reason=="included",ID])) stop("Anchor failed frequency check; do not relax thresholds automatically")
  ma <- z[reason=="included",.(SNP=ID,A1,A2,freq=FREQ,b=BETA,se=SE,p=P,N=NTOT)]
  # P=0 can be finite-number underflow. Floor only for the GCTA text interface.
  ma[, p := pmax(p,.Machine$double.xmin)]
  mafile <- file.path(od,"PTSD_NTOT_COJO_only.ma")
  data.table::fwrite(ma,mafile,sep="\t")
  dput(list(build="GRCh37",region=c(inp$lower,inp$upper),reference_n=503,
       diff_freq=.10,cojo_collinear=.90,N="source NTOT",outcome="case-control PTSD log OR",
       gcta=gcta,interpretation="exploratory conditional SNP association"),file=file.path(od,"parameters.txt"))
  runs <- data.table::data.table(conditioned_gene=c("FES","FURIN"),
      conditioned_snp=c("rs1894401","rs2071410"),target_gene=c("FURIN","FES"),
      target_snp=c("rs2071410","rs1894401"))
  result <- vector("list",nrow(runs))
  for(i in seq_len(nrow(runs))) {
    rr <- runs[i]; snplist <- file.path(od,paste0("condition_",rr$conditioned_gene,".snplist"))
    writeLines(rr$conditioned_snp,snplist)
    pref <- file.path(od,paste0("ptsd_cond_",rr$conditioned_gene))
    ptsd_followup_note("COJO: PTSD conditioned on ",rr$conditioned_gene)
    err <- tryCatch({
      ptsd_run_external(gcta,c("--bfile",refout,"--chr","15","--cojo-file",mafile,
          "--cojo-cond",snplist,"--diff-freq","0.10","--cojo-collinear","0.90",
          "--out",pref,"--thread-num","2"),paste0(pref,"_console.txt")); NULL
    },error=function(e) conditionMessage(e))
    u <- ma[SNP==rr$target_snp]
    marginal <- data.table::data.table(outcome_A1=u$A1,outcome_A2=u$A2,
          beta_marginal=u$b,se_marginal=u$se,p_marginal=u$p,N=u$N)
    new <- if(is.null(err)) ptsd_cojo_effect(paste0(pref,".cma.cojo"),rr$target_snp,u$A1,u$A2)
           else data.table::data.table(status=paste("program_error",err))
    oldfile <- file.path(P$project_root,"results/cojo_chr15",paste0("ptsd_cond_",rr$conditioned_gene,".cma.cojo"))
    oldres <- ptsd_cojo_effect(oldfile,rr$target_snp,u$A1,u$A2)
    data.table::setnames(oldres,names(oldres),paste0("old_",names(oldres)))
    result[[i]] <- cbind(rr,marginal,new,oldres)
  }
  ans <- data.table::rbindlist(result,fill=TRUE)
  ans[, interpretation := "exploratory conditional SNP association; no protein causal estimate"]
  data.table::fwrite(ans,file.path(od,"COJO_PTSD_comparison.csv"))
  print(ans)
  if(any(ans$status!="estimated")) stop("Some COJO target effects were not estimable; see comparison and raw logs")
  writeLines("COMPUTATION COMPLETE; small-reference and scientific review still required.",file.path(od,"SUCCESS.txt"))
  invisible(od)
}

run_ptsd_sem_neff <- function(P) {
  root <- P$project_root
  out <- file.path(P$followup_out, "SEM_NEFF")
  if (!dir.exists(file.path(root, "results"))) stop("Missing D:/PTSD/results")
  if (!dir.create(out)) stop("Could not create new output folder")
  oldwd <- getwd()
  out_depth <- sink.number()
  msg_depth <- sink.number(type = "message")
  con <- file(file.path(out, "RUN.log"), open = "wt")
  sink(con, split = TRUE)
  sink(con, type = "message")
  on.exit({
    if (sink.number(type = "message") != msg_depth) sink(type = "message")
    while (sink.number() > out_depth) sink()
    close(con)
    setwd(oldwd)
    cat("\nOutput directory: ", out, "\n", sep = "")
  }, add = TRUE)
  tryCatch({
    cat("Start: ", format(Sys.time()), "\nOutput: ", out, "\n", sep = "")
    pkgs <- c("GenomicSEM", "data.table", "R.utils")
    missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing)) stop("Missing installed packages: ", paste(missing, collapse = ", "))
    suppressPackageStartupMessages(library(data.table))
    writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_start.txt"))
    writeLines(capture.output(packageDescription("GenomicSEM")), file.path(out, "GenomicSEM_description.txt"))
    writeLines(capture.output(dput(run_ptsd_sem_neff)), file.path(out, "rerun_function_snapshot.R.txt"))
    if (!"Q_SNP" %in% names(formals(GenomicSEM::userGWAS))) stop("Installed userGWAS does not expose Q_SNP. Send package description for review.")
    # Preserve the installed implementations, including development commit details.
    for (fn in c("munge", "ldsc", "sumstats", "usermodel", "userGWAS")) {
      writeLines(capture.output(dput(getExportedValue("GenomicSEM", fn))),
                 file.path(out, paste0("installed_", fn, ".R.txt")))
    }

    cfg <- new.env(parent = globalenv()); cfg$P <- P
    config_file <- file.path(P$package_root, "00_config.R")
    if (file.exists(config_file)) sys.source(config_file, envir = cfg)
    pick <- function(key, fallback) {
      if (exists(key, envir = cfg, inherits = FALSE)) {
        x <- get(key, envir = cfg, inherits = FALSE)
        if (is.character(x) && length(x) == 1L && !grepl("/path/to/", x, fixed = TRUE)) return(x)
      }
      fallback
    }
    inputs <- c(
      PTSD = P$ptsd_sem_ma,
      MDD = pick("mdd_path", file.path(root, "GWAS/PGC_MDD_2025/27061255/daner/daner_pgc_mdd_no23andMe-noUKBB_eur_hg19_v3.49.24.11.neff.gz")),
      ANX = pick("anx_path", file.path(root, "GWAS/PGC_Anxiety_2026/31389910/ANX_2026_daner_fullANX_v12_woUTAH_11022026.gz")),
      pQTL = pick("pqtl_path", file.path(root, "pQTL/ukbppp_cis_pqtl.rds")),
      HM3 = pick("hm3_path", P$hm3),
      reference = pick("ref_sumstats", P$sumstats_reference)
    )
    ld_dir <- pick("ldsc_dir", P$ldsc)
    ld_dir <- paste0(sub("[/\\\\]+$", "", ld_dir), "/")
    missing_files <- inputs[!file.exists(inputs)]
    if (length(missing_files)) stop("Input file(s) missing:\n", paste(names(missing_files), missing_files, collapse = "\n"))
    ld_files <- unlist(lapply(1:22, function(i) paste0(ld_dir, i, c(".l2.ldscore.gz", ".l2.M_5_50"))))
    if (any(!file.exists(ld_files))) stop("LD reference missing:\n", paste(ld_files[!file.exists(ld_files)], collapse = "\n"))
    # Enforce the intended phenotype releases; a different file needs review.
    if (basename(inputs["MDD"]) != "daner_pgc_mdd_no23andMe-noUKBB_eur_hg19_v3.49.24.11.neff.gz") stop("MDD input differs from the intended no23andMe-noUKBB release.")
    if (basename(inputs["ANX"]) != "ANX_2026_daner_fullANX_v12_woUTAH_11022026.gz") stop("ANX input differs from the intended woUTAH release.")
    setwd(out)
    manifest_paths <- c(unname(inputs), ld_files, if (file.exists(config_file)) config_file)
    cat("Hashing input files...\n")
    manifest <- data.frame(path = normalizePath(manifest_paths, winslash = "/"),
                           bytes = file.info(manifest_paths)$size,
                           md5 = unname(tools::md5sum(manifest_paths)))
    write.csv(manifest, "input_manifest.csv", row.names = FALSE)
    traits <- c("PTSD", "MDD", "ANX")
    sample_prev <- c(0.5, 0.5, 0.5)
    population_prev <- c(0.07, 0.15, 0.10)
    n_vec <- c(NA_real_, 829250, 390020)
    pars <- list(trait_order = traits, N = n_vec, sample.prev = sample_prev,
                 population.prev = population_prev, munge.info = 0.9,
                 sumstats.info = 0.6, MAF = 0.01, estimation = "DWLS",
                 note = "SENSITIVITY ONLY: PTSD uses the source per-SNP NEFF column with sample prevalence 0.5. Separate count-pattern diagnostics support a cohort-additive structure, without reconstructing all absolute cohort weights (S24c). MDD/ANX settings and the factor model are retained from the reviewed rerun; this is not a new audit of their sample-size definitions.")
    dput(pars, file = "parameters.txt")
    print(pars)
    for (nm in c("MDD", "ANX")) {
      hdr <- fread(inputs[nm], nrows = 5)
      writeLines(capture.output(print(hdr)), paste0(nm, "_input_preview.txt"))
    }

    cat("\nSTEP 1/5: Prepare and munge source summary statistics\n")
    ptsd <- fread(inputs["PTSD"])
    need <- c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")
    if (!all(need %in% names(ptsd))) stop("PTSD .ma requires: ", paste(need, collapse = ", "))
    if (anyDuplicated(ptsd$SNP)) stop("Duplicate PTSD SNP IDs: resolve before running.")
    if (any(!is.finite(ptsd$N) | ptsd$N <= 0)) stop("Invalid PTSD per-SNP N")
    writeLines(capture.output(summary(ptsd$N)), "PTSD_N_summary.txt")
    setnames(ptsd, need, c("SNP", "A1", "A2", "MAF", "BETA", "SE", "P", "N"))
    prep <- file.path(out, "PTSD_prepared.txt")
    fwrite(ptsd, prep, sep = "\t")
    rm(ptsd); gc()
    gwas_files <- c(prep, unname(inputs[c("MDD", "ANX")]))
    GenomicSEM::munge(files = gwas_files, hm3 = unname(inputs["HM3"]),
                     trait.names = traits, N = n_vec, info.filter = 0.9, maf.filter = 0.01)
    munged <- paste0(traits, ".sumstats.gz")
    if (!all(file.exists(munged))) stop("Munge did not create all three files.")
    write.csv(data.frame(file = munged, md5 = unname(tools::md5sum(munged))),
              "munged_manifest.csv", row.names = FALSE)

    # Verify what munge actually used, not only the requested setting.
    n_source <- fread(inputs["PTSD"], select = c("SNP", "N"))
    n_munged <- fread("PTSD.sumstats.gz", select = c("SNP", "N"))
    n_check <- merge(n_munged, n_source, by = "SNP", all.x = TRUE,
                     suffixes = c("_munged", "_source"))
    n_bad <- n_check[!is.finite(N_source) | !is.finite(N_munged) |
                     abs(N_source - N_munged) > 1e-6]
    fwrite(n_bad, "PTSD_munged_N_mismatches.csv")
    if (nrow(n_bad)) stop("Munged PTSD N does not match source NEFF")
    n_audit <- rbindlist(lapply(seq_along(traits), function(j) {
      v <- fread(munged[j], select = "N")$N
      data.table(trait = traits[j], n_snps = length(v), min_N = min(v),
                 median_N = median(v), max_N = max(v), sample_prevalence = sample_prev[j])
    }))
    fwrite(n_audit, "MUNGE_N_audit.csv")
    rm(n_source, n_munged, n_check, n_bad); gc()

    cat("\nSTEP 2/5: Fresh LDSC and trait models\n")
    ld <- GenomicSEM::ldsc(traits = munged, sample.prev = sample_prev,
                          population.prev = population_prev, ld = ld_dir, wld = ld_dir,
                          trait.names = traits)
    saveRDS(ld, "LDSCoutput_stage0.rds")
    S <- as.matrix(ld$S)
    write.csv(S, "LDSC_S.csv")
    write.csv(as.matrix(ld$V), "LDSC_V.csv")
    write.csv(cov2cor(S), "genetic_correlations.csv")
    if (!all(is.finite(S)) || any(diag(S) <= 0)) stop("Invalid LDSC covariance matrix.")
    eig <- eigen((S+t(S))/2, symmetric=TRUE, only.values=TRUE)$values
    writeLines(capture.output(eig), "LDSC_S_eigenvalues.txt")
    if (min(eig) <= 0) warning("LDSC S is not positive definite: model results need explicit review.")
    # Compare the package's own liability-scale log against the returned matrix.
    logs <- list.files(out, pattern = "_ldsc\\.log$", full.names = TRUE)
    if (length(logs) == 1L) {
      lines <- readLines(logs, warn = FALSE)
      hl <- grep("Total Liability Scale h2:", lines, value = TRUE, fixed = TRUE)
      if (length(hl) == 3L) {
        h <- as.numeric(sub(".*Total Liability Scale h2: *([0-9.eE+-]+).*", "\\1", hl))
        check <- data.frame(trait = traits, log_h2 = h, matrix_h2 = diag(S), difference = diag(S) - h)
        write.csv(check, "LDSC_log_matrix_check.csv", row.names = FALSE)
        if (any(!is.finite(h)) || any(abs(h - diag(S)) > 0.00006)) stop("Fresh LDSC log and S matrix disagree; inspect LDSC_log_matrix_check.csv before continuing.")
      } else warning("Could not parse three liability h2 values; check the LDSC log manually.")
    } else warning("Could not identify a unique LDSC log; check manually.")
    model_free <- 'F1 =~ NA*PTSD + MDD + ANX
                   F1 ~~ 1*F1
                   PTSD ~~ PTSD
                   MDD ~~ MDD
                   ANX ~~ ANX'
    model_fixed <- 'F1 =~ NA*PTSD + MDD + ANX
                    F1 ~~ 1*F1
                    PTSD ~~ PTSD
                    MDD ~~ 0*MDD
                    ANX ~~ ANX'
    writeLines(model_fixed, "model_primary.txt")
    free_fit <- tryCatch(GenomicSEM::usermodel(covstruc = ld, estimation = "DWLS", model = model_free,
                         CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE),
                         error = function(e) list(error = conditionMessage(e)))
    fit <- GenomicSEM::usermodel(covstruc = ld, estimation = "DWLS", model = model_fixed,
                                    CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
    saveRDS(list(free = free_fit, constrained = fit), "stage0b_models.rds")
    writeLines(capture.output(print(fit)), "constrained_model_full.txt")
    writeLines(capture.output(print(free_fit)), "unconstrained_model_full.txt")
    if (is.data.frame(fit$results)) fwrite(fit$results, "constrained_model_parameters.csv")
    if (is.data.frame(fit$modelfit)) fwrite(fit$modelfit, "constrained_model_fit.csv")

    cat("\nSTEP 3/5: Fresh sumstats harmonization (may be slow)\n")
    ss <- GenomicSEM::sumstats(files = gwas_files, ref = unname(inputs["reference"]),
              trait.names = traits, se.logit = rep(TRUE, 3), OLS = rep(FALSE, 3),
              linprob = rep(FALSE, 3), N = n_vec, info.filter = 0.6, maf.filter = 0.01)
    saveRDS(ss, "stage1_sumstats.rds")
    ss <- as.data.table(ss)
    if (!all(c("SNP", "A1", "A2") %in% names(ss)) || anyDuplicated(ss$SNP)) stop("Invalid/duplicate harmonized SNP IDs")
    pqtl <- as.data.table(readRDS(inputs["pQTL"]))
    cols <- c("protein", "SNP", "beta", "se", "effect_allele", "other_allele")
    if (!all(cols %in% names(pqtl))) stop("Missing required pQTL columns: ", paste(setdiff(cols, names(pqtl)), collapse = ", "))
    if (anyDuplicated(pqtl$protein)) stop("More than one instrument per protein: resolve selection first.")
    ss_inst <- ss[SNP %in% pqtl$SNP]
    rm(ss); gc()
    if (!nrow(ss_inst)) stop("No direct instruments in fresh sumstats")
    fwrite(ss_inst, "instrument_sumstats.csv")
    cat("Direct instrument SNPs: ", nrow(ss_inst), "\n", sep = "")

    cat("\nSTEP 4/5: Fresh instrument-only userGWAS; no cache\n")
    if (!"Q_SNP" %in% names(formals(GenomicSEM::userGWAS))) stop("Installed userGWAS does not expose Q_SNP. Send package description for review.")
    gw <- GenomicSEM::userGWAS(covstruc = ld, SNPs = as.data.frame(ss_inst),
          estimation = "DWLS", model = paste(model_fixed, "F1 ~ SNP", sep = "\n"),
          sub = "F1~SNP", Q_SNP = TRUE, cores = 1, parallel = FALSE, printwarn = TRUE)
    saveRDS(gw, "stage12_factor_snpeffects_userGWAS.rds")
    fg <- as.data.table(gw[[1]])
    fwrite(fg, "userGWAS_raw.csv")
    required <- c("SNP", "est", "SE", "Q_SNP", "Q_SNP_pval")
    if (!all(required %in% names(fg))) stop("Unexpected userGWAS columns; inspect userGWAS_raw.csv")
    if (anyDuplicated(fg$SNP)) stop("Duplicate userGWAS SNP rows")
    if ("error" %in% names(fg)) {
      bad <- !is.na(fg$error) & !as.character(fg$error) %in% c("0", "FALSE", "")
      if (any(bad)) stop("userGWAS contains model errors; inspect raw output")
    }
    if (any(!is.finite(fg$est) | !is.finite(fg$SE) | fg$SE <= 0 |
            !is.finite(fg$Q_SNP_pval) | fg$Q_SNP_pval < 0 | fg$Q_SNP_pval > 1)) stop("Invalid userGWAS estimates or Q P values")

    cat("\nSTEP 5/5: Full allele-pair matching, Wald MR and BH\n")
    key <- fg[, .(SNP, factor_beta_raw = est, se_factor = SE, Q = Q_SNP, Q_pval = Q_SNP_pval)]
    key <- merge(key, ss_inst[, .(SNP, ssA1 = toupper(A1), ssA2 = toupper(A2))], by = "SNP")
    a <- merge(pqtl, key, by = "SNP", all.x = TRUE)
    a[, `:=`(EA = toupper(effect_allele), OA = toupper(other_allele))]
    a[, direct := EA == ssA1 & OA == ssA2]
    a[, reverse := EA == ssA2 & OA == ssA1]
    a[, palindromic := paste0(EA, OA) %in% c("AT", "TA", "CG", "GC")]
    a[, reason := "included"]
    a[is.na(factor_beta_raw), reason := "No direct instrument in factor output"]
    a[reason == "included" & !(direct %in% TRUE | reverse %in% TRUE), reason := "Allele pair mismatch"]
    a[reason == "included" & (!is.finite(beta) | beta == 0 | !is.finite(se) | se <= 0), reason := "Invalid pQTL effect or SE"]
    # Palindromes are not resolved from allele labels alone. Exclude and report.
    a[reason == "included" & palindromic, reason := "Palindromic SNP: requires EAF/strand review"]
    fwrite(a, "instrument_audit.csv")
    fwrite(a[reason != "included", .(protein, SNP, reason)], "stage12_not_testable.csv")
    r <- a[reason == "included"]
    if (!nrow(r)) stop("No valid harmonized protein instruments")
    r[, beta_factor := factor_beta_raw * ifelse(direct, 1, -1)]
    r[, mr_beta := beta_factor / beta]
    r[, mr_se := abs(se_factor / beta)]
    r[, mr_z := mr_beta / mr_se]
    r[, mr_p := 2 * pnorm(-abs(mr_z))]
    r[, `:=`(mr_fdr = p.adjust(mr_p, "BH"), Q_fdr = p.adjust(Q_pval, "BH"))]
    if (!"cis_gene" %in% names(r)) r[, cis_gene := NA_character_]
    ans <- r[, .(protein, cis_gene, outcome_snp = SNP,
          exposure_EA = EA, exposure_OA = OA, sumstats_A1 = ssA1, sumstats_A2 = ssA2,
          pqtl_beta = beta, pqtl_se = se, beta_factor, se_factor,
          mr_beta, mr_se, mr_z, mr_p, mr_fdr, Q, Q_pval, Q_fdr)]
    setorder(ans, mr_p)
    fwrite(ans, "stage12_factorMR_results.csv")
    candidates <- c("AKT3", "CD40", "CGREF1", "FES", "FURIN", "SIRPA", "CD101", "KHK", "SNX18", "UBE2L6")
    fwrite(ans[protein %in% candidates], "candidate_results.csv")
    counts <- data.frame(testable_proteins = nrow(ans), unique_SNPs = uniqueN(ans$outcome_snp),
                        MR_FDR_005 = sum(ans$mr_fdr < .05), Q_nominal_005 = sum(ans$Q_pval < .05),
                        Q_FDR_005 = sum(ans$Q_fdr < .05))
    write.csv(counts, "summary_counts.csv", row.names = FALSE)
    print(counts)
    print(ans[protein %in% candidates])
    legacy_path <- file.path(P$package_root, "results_snapshot/factor_MR.csv")
    if (file.exists(legacy_path)) {
      legacy <- fread(legacy_path)
      cmp <- merge(ans[, .(protein, outcome_snp, exposure_EA, exposure_OA, mr_beta, mr_se, mr_p, mr_fdr, Q_pval, Q_fdr)],
                   legacy[, .(protein, outcome_snp, exposure_EA, exposure_OA, mr_beta, mr_se, mr_p, mr_fdr, Q_pval, Q_fdr)],
                   by = c("protein", "outcome_snp"), all = TRUE, suffixes = c("_NEFF", "_reviewed"))
      cmp[, allele_agreement := exposure_EA_NEFF == exposure_EA_reviewed &
                                  exposure_OA_NEFF == exposure_OA_reviewed]
      cmp[, sign_changed := ifelse(allele_agreement %in% TRUE,
                                  sign(mr_beta_NEFF) != sign(mr_beta_reviewed), NA)]
      cmp[, MR_FDR_status_changed := (mr_fdr_NEFF < .05) != (mr_fdr_reviewed < .05)]
      cmp[, Q_FDR_status_changed := (Q_fdr_NEFF < .05) != (Q_fdr_reviewed < .05)]
      fwrite(cmp, "comparison_with_reviewed.csv")
      fwrite(cmp[protein %in% candidates], "candidate_comparison.csv")
      old_counts <- data.table(analysis="reviewed_20261009", testable_proteins=nrow(legacy),
          unique_SNPs=uniqueN(legacy$outcome_snp), MR_FDR_005=sum(legacy$mr_fdr < .05),
          Q_nominal_005=sum(legacy$Q_pval < .05), Q_FDR_005=sum(legacy$Q_fdr < .05))
      new_counts <- as.data.table(counts); new_counts[, analysis := "NEFF_sensitivity"]
      fwrite(rbindlist(list(old_counts,new_counts), use.names=TRUE), "counts_comparison.csv")
    } else {
      warning("Reviewed factor_MR.csv is missing; new results saved without an automatic comparison.")
    }
    # Verify source inputs were not modified during the run.
    if (!identical(unname(tools::md5sum(manifest_paths)), manifest$md5)) stop("Input files changed during run; do not use results.")
    writeLines(capture.output(sessionInfo()), "sessionInfo_end.txt")
    writeLines(c("COMPUTATION COMPLETE; NEFF-definition and scientific review still required. Not a submission-final replacement.", format(Sys.time())), "SUCCESS.txt")
    cat("\nCOMPUTATION COMPLETE. Upload summary_counts.csv, candidate_results.csv,\n",
        "candidate_comparison.csv, MUNGE_N_audit.csv, constrained_model_full.txt and RUN.log for review.\n")
  }, error = function(e) {
    writeLines(c(format(Sys.time()), conditionMessage(e)), file.path(out, "FAILED.txt"))
    cat("\nSTOPPED: ", conditionMessage(e), "\n", sep = "")
    stop(e)
  })
  invisible(out)
}

ptsd_followup_20261016 <- function() {
  if (!requireNamespace("data.table",quietly=TRUE)) stop("Missing package: data.table")
  if (!requireNamespace("R.utils",quietly=TRUE)) stop("Missing package: R.utils")
  suppressPackageStartupMessages(library(data.table))
  package_root <- getOption("ptsd.followup.package_root",getOption("ptsd.current.root"))
  paths <- file.path(package_root,"00_paths.R")
  if(!file.exists(paths)) stop("Current package is missing: ",paths)
  old_options <- list(ptsd.current.root=getOption("ptsd.current.root"),ptsd.current.paths=getOption("ptsd.current.paths"))
  options(ptsd.current.root=package_root)
  on.exit(options(old_options),add=TRUE)
  env <- new.env(parent=globalenv()); sys.source(paths,envir=env)
  P <- env$P; P$package_root <- package_root
  P$mode <- getOption("ptsd.followup.mode","both")
  if(length(P$mode)!=1L || !P$mode%in%c("both","sem","cojo")) stop("mode must be both, sem, or cojo")
  # P$source_casecontrol is defined by this bundle's 00_paths.R.
  dir.create(P$result_dir,recursive=TRUE,showWarnings=FALSE)
  P$followup_out <- tempfile(paste0("sem_cojo_followup_",format(Sys.time(),"%Y%m%d_%H%M%S"),"_"),tmpdir=P$result_dir)
  if(!dir.create(P$followup_out)) stop("Cannot create new output directory")
  review <- file.path(P$followup_out,"SEND_BACK"); dir.create(review)
  ptsd_followup_note("Output: ",P$followup_out)
  writeLines("2026-10-17: PLINK2 force-intersect and complete SEND_BACK logs",
             file.path(P$followup_out,"SCRIPT_REVISION.txt"))
  writeLines(capture.output(sessionInfo()),file.path(P$followup_out,"sessionInfo.txt"))
  writeLines(capture.output(dput(P)),file.path(P$followup_out,"paths_and_mode.txt"))
  status <- data.table(module=character(),status=character(),detail=character())
  add_status <- function(module,state,detail) {
    status <<- rbind(status,data.table(module=module,status=state,detail=detail))
    fwrite(status,file.path(P$followup_out,"FOLLOWUP_STATUS.csv"))
  }
  on.exit({
    # Collect only small review outputs, never GWAS data or genotype files.
    selected <- c("FOLLOWUP_STATUS.csv","INTERPRETATION_LIMITS.txt","sessionInfo.txt","SCRIPT_REVISION.txt",
       "INPUTS/PTSD_sample_size_summary.csv","INPUTS/NEFF_diagnostics.csv","INPUTS/candidate_source_rows.csv",
       "INPUTS/PTSD_full_source_comparison.csv",
       "COJO_PTSD/LD_reference_audit.csv","COJO_PTSD/GCTA_candidates.csv",
       "COJO_PTSD/reference_anchor_rows.csv","COJO_PTSD/anchor_frequency_check.csv",
       "COJO_PTSD/COJO_PTSD_comparison.csv","COJO_PTSD/FAILED.txt",
       "COJO_PTSD/parameters.txt",
       "COJO_PTSD/PLINK_extract_console.txt","COJO_PTSD/PLINK_extract_console.txt.command.txt",
       "COJO_PTSD/EUR503_chr15.log",
       "COJO_PTSD/PLINK_frequency_console.txt","COJO_PTSD/PLINK_frequency_console.txt.command.txt",
       "COJO_PTSD/EUR503_frequency.log",
       "COJO_PTSD/ptsd_cond_FES.log","COJO_PTSD/ptsd_cond_FURIN.log",
       "COJO_PTSD/ptsd_cond_FES_console.txt","COJO_PTSD/ptsd_cond_FURIN_console.txt",
       "COJO_PTSD/ptsd_cond_FES_console.txt.command.txt","COJO_PTSD/ptsd_cond_FURIN_console.txt.command.txt",
       "SEM_NEFF/summary_counts.csv","SEM_NEFF/candidate_results.csv","SEM_NEFF/counts_comparison.csv",
       "SEM_NEFF/stage12_factorMR_results.csv","SEM_NEFF/comparison_with_reviewed.csv",
       "SEM_NEFF/candidate_comparison.csv","SEM_NEFF/MUNGE_N_audit.csv","SEM_NEFF/LDSC_log_matrix_check.csv",
       "SEM_NEFF/LDSC_S.csv","SEM_NEFF/genetic_correlations.csv","SEM_NEFF/constrained_model_full.txt",
       "SEM_NEFF/unconstrained_model_full.txt","SEM_NEFF/parameters.txt","SEM_NEFF/RUN.log","SEM_NEFF/FAILED.txt")
    for(rel in selected) {
      src <- file.path(P$followup_out,rel)
      if(file.exists(src)) file.copy(src,file.path(review,gsub("/","__",rel,fixed=TRUE)),overwrite=FALSE)
    }
    cat("\n========================================\n")
    print(status)
    cat("\nNEW OUTPUT: ",P$followup_out,"\nSEND BACK the small files in: ",review,"\n",sep="")
    cat("Do not substitute these results into the manuscript or frozen reports until reviewed.\n")
  },add=TRUE)
  inp <- tryCatch(ptsd_prepare_followup(P),error=function(e) {
    add_status("source_preparation","FAILED",conditionMessage(e)); NULL
  })
  if(is.null(inp)) return(invisible(P$followup_out))
  P$ptsd_sem_ma <- inp$sem_ma
  add_status("source_preparation","COMPLETE","Official case-control source verified; diagnostic outputs saved")
  if(P$mode%in%c("both","cojo")) {
    tryCatch({
      run_ptsd_cojo_followup(P,inp)
      add_status("COJO_PTSD","COMPLETE_REVIEW_REQUIRED","Small-reference sensitivity; trauma COJO remains pending")
    },error=function(e) {
      dir.create(file.path(P$followup_out,"COJO_PTSD"),showWarnings=FALSE)
      writeLines(conditionMessage(e),file.path(P$followup_out,"COJO_PTSD/FAILED.txt"))
      add_status("COJO_PTSD","FAILED_OR_BLOCKED",conditionMessage(e))
      ptsd_followup_note("COJO stopped: ",conditionMessage(e),"; proceeding to any requested SEM step.")
    })
  }
  rm(inp); gc()
  if(P$mode%in%c("both","sem")) {
    tryCatch({
      run_ptsd_sem_neff(P)
      add_status("SEM_NEFF","COMPLETE_REVIEW_REQUIRED","Source NEFF + sample prevalence 0.5; separate S24c audit supports its information scale, with incomplete absolute cohort weights")
    },error=function(e) {
      add_status("SEM_NEFF","FAILED_OR_BLOCKED",conditionMessage(e))
      ptsd_followup_note("SEM stopped: ",conditionMessage(e))
    })
  }
  writeLines(c("Computations attempted; see FOLLOWUP_STATUS.csv for actual per-module success.",
               "Outputs are sensitivity results awaiting review, not a final script/manuscript release."),
             file.path(P$followup_out,"READ_ME_FIRST.txt"))
  invisible(P$followup_out)
}

if (!isTRUE(getOption("ptsd.followup.load_only",FALSE))) ptsd_followup_20261016()
