# Current package: path configuration only; scientific calculations preserved.
P <- getOption("ptsd.current.paths")
if (is.null(P)) stop("First source 00_START.R; then use ptsd_run().")
# FES/FURIN regional colocalization audit, 2026-10-08
# Source this script after extracting the entire ZIP into one directory.
# Requires data.table, R.utils, coloc, susieR; original Windows inputs are read-only.
# Exploratory reanalysis: external LD, variable case-control N, no automatic causal claim.
(function() {
  sf <- file.path(P$package_root,"scripts","PTSD_FES_FURIN_coloc_rerun.R")
  here <- P$package_root
  plink <- P$plink2
  ld <- P$ld_prefix
  base <- P$protein_raw_base
  paths <- c(FES=file.path(base,"FES/discovery_chr15_FES_P07332_OID21207_v1_Oncology.gz"),
             FURIN=file.path(base,"FURIN/discovery_chr15_FURIN_P09958_OID21514_v1_Oncology.gz"))
  ypath <- P$ptsd_region
  required <- c("data.table","R.utils","coloc","susieR")
  absent <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
  if(length(absent)) stop("Install packages first: install.packages(c(",paste(shQuote(absent),collapse=","),"))")
  library(data.table)
  inputs <- c(plink,paste0(ld,c(".bed",".bim",".fam")),paths,ypath)
  if(any(!file.exists(inputs))) stop("Missing: ",paste(inputs[!file.exists(inputs)],collapse="; "))
  out <- tempfile(paste0("coloc_FES_FURIN_",format(Sys.time(),"%Y%m%d_%H%M%S"),"_"),tmpdir=P$result_dir)
  dir.create(out,recursive=TRUE)
  cat("Output directory: ",out,"\n",sep="")
  log <- file(file.path(out,"RUN.log"),"wt")
  sink(log,split=TRUE); sink(log,type="message")
  on.exit({sink(type="message");sink();close(log)},add=TRUE)
  file.copy(sf,file.path(out,"script_snapshot.R.txt"))
  writeLines(c("Region: each sentinel +/-500kb in GRCh37; position from protein ID, not GENPOS.",
    "No association P-value filter. Biallelic A/C/G/T only; protein INFO>=0.8; reference MAF>=0.01; genotype missing<=0.05.",
    "Exact allele-pair matching in either order; no automatic strand-complement rescue.",
    "All effects aligned to BIM B1, explicitly counted by PLINK --export-allele.",
    "Frequency QC: aligned exposure/outcome vs reference difference<=0.15; palindromes also require direct frequency distance<=0.10 and opposite orientation worse by>0.10.",
    "LD: Pearson correlation of reference genotypes after per-variant mean imputation; no clumping; no automatic ridge or SNP removal based on z/LD diagnostics.",
    "SuSiE: L=10, coverage=.95, min_abs_corr=.5, max_iter=1000; residual variance fixed for external LD.",
    "Protein N: observed discovery sample size. PTSD primary approximation: median NTOT of retained SNPs; sensitivity: median reported NEFF.",
    "Scalar N approximates variable-N case-control summary statistics; neither choice resolves cohort overlap or ascertainment. Review both fits and diagnostics.",
    "ABF: single-causal-variant sensitivity analysis; protein sdY estimated by coloc using beta/SE, MAF and N; outcome type cc with log-odds beta/SE.",
    "Prior grid p12=1e-7,5e-7,1e-6,5e-6,1e-5,5e-5,1e-4; p1=p2=1e-4.",
    "Missing credible sets, nonconvergence or LD inconsistency are not evidence against colocalization.",
    "PTSD subset is restricted to LD-reference SNPs; this does not establish complete original GWAS regional coverage."),file.path(out,"parameters.txt"))
  tryCatch({
    cat("STEP 1/4: Inputs and regional allele matching\n")
    # Hash summary inputs; large reference BED metadata retained without hashing.
    fi <- file.info(inputs)
    fwrite(data.table(path=inputs,bytes=fi$size,mtime=as.character(fi$mtime)),file.path(out,"input_manifest.csv"))
    print(tools::md5sum(c(paths,ypath)))
    print(system2(plink,"--version",stdout=TRUE,stderr=TRUE))
    b <- fread(paste0(ld,".bim"),header=FALSE,col.names=c("CHR","SNP","CM","BP","B1","B2"))
    b <- b[CHR==15 & BP>=90920940 & BP<=91929042]
    b <- b[!duplicated(SNP) & !duplicated(SNP,fromLast=TRUE)]
    y <- fread(ypath)
    stopifnot(!anyDuplicated(y$SNP),all(y$NTOT==y$NCASE+y$NCON))
    setnames(y,c("A1","A2","FREQ","BETA","SE","P","BP"),c("YEA","YOA","YFREQ","YBETA","YSE","YP","YBP"))
    all_data <- list(); counts <- list(); statuses <- list(); abfs <- list(); susies <- list()
    centers <- c(FES=91429042,FURIN=91420940)
    sentinels <- c(FES="rs1894401",FURIN="rs2071410")
    for(g in names(paths)) {
      p <- fread(paths[[g]])
      parts <- tstrsplit(p$ID,":",fixed=TRUE)
      p[, BP:=as.integer(parts[[2]])]
      p <- p[CHROM==15 & abs(BP-centers[[g]])<=500000 & TEST=="ADD"]
      p[, pair:=paste(pmin(ALLELE0,ALLELE1),pmax(ALLELE0,ALLELE1),sep="/")]
      bb <- copy(b); bb[,pair:=paste(pmin(B1,B2),pmax(B1,B2),sep="/")]
      d <- merge(p,bb,by=c("BP","pair"),allow.cartesian=TRUE)
      d <- merge(d,y,by="SNP",allow.cartesian=TRUE)
      d[, reason:="keep"]
      d[duplicated(SNP)|duplicated(SNP,fromLast=TRUE),reason:="duplicate_mapping"]
      d[reason=="keep" & !(ALLELE0 %in% c("A","C","G","T") & ALLELE1 %in% c("A","C","G","T") & ALLELE0!=ALLELE1),reason:="not_biallelic_SNV"]
      d[reason=="keep" & !(BP==YBP & ((YEA==B1 & YOA==B2)|(YEA==B2 & YOA==B1))),reason:="outcome_pair_or_position_mismatch"]
      d[reason=="keep" & !(is.finite(BETA)&is.finite(SE)&SE>0&is.finite(YBETA)&is.finite(YSE)&YSE>0&is.finite(INFO)&INFO>=0.8&N>0&NTOT>0&NEFF>0&A1FREQ>0&A1FREQ<1&YFREQ>0&YFREQ<1),reason:="invalid_statistics_or_INFO"]
      d[, `:=`(xb=ifelse(ALLELE1==B1,BETA,-BETA),yb=ifelse(YEA==B1,YBETA,-YBETA),
                xf=ifelse(ALLELE1==B1,A1FREQ,1-A1FREQ),yf=ifelse(YEA==B1,YFREQ,1-YFREQ),
                palindrome=pair %in% c("A/T","C/G"))]
      fwrite(d,file.path(out,paste0(g,"_preLD_audit.csv")))
      counts[[g]] <- data.table(protein=g,raw_region=nrow(p),matched_before_QC=nrow(d),preLD_pass=sum(d$reason=="keep"))
      all_data[[g]] <- d[reason=="keep"]
    }
    ids <- unique(rbindlist(lapply(all_data,function(d)d[,.(SNP,B1)])))
    stopifnot(!anyDuplicated(ids$SNP),nrow(ids)>10)
    fwrite(ids[,.(SNP)],file.path(out,"extract.txt"),col.names=FALSE)
    fwrite(ids,file.path(out,"counted_alleles.txt"),sep="\t",col.names=FALSE)
    cat("STEP 2/4: Export reference genotypes\n")
    pref <- file.path(out,"reference")
    args <- c("--bfile",shQuote(ld),"--extract",shQuote(file.path(out,"extract.txt")),
              "--export","A","--export-allele",shQuote(file.path(out,"counted_alleles.txt")),"--out",shQuote(pref))
    writeLines(paste(shQuote(plink),paste(args,collapse=" ")),file.path(out,"plink_command.txt"))
    code <- system2(plink,args,stdout=file.path(out,"plink_console.txt"),stderr=file.path(out,"plink_stderr.txt"))
    if(code!=0 || !file.exists(paste0(pref,".raw"))) stop("PLINK genotype export failed; see plink logs")
    gr <- fread(paste0(pref,".raw"))
    gn <- paste(ids$SNP,ids$B1,sep="_")
    if(!all(gn %in% names(gr))) stop("Counted allele headers do not match requested alleles")
    G <- as.matrix(gr[,..gn]); storage.mode(G)<-"double"; colnames(G)<-ids$SNP
    miss <- colMeans(is.na(G)); fr <- colMeans(G,na.rm=TRUE)/2
    for(j in seq_len(ncol(G))) G[is.na(G[,j]),j] <- 2*fr[j]
    good <- is.finite(fr)&pmin(fr,1-fr)>=.01&miss<=.05
    G <- G[,good,drop=FALSE]
    ref <- data.table(SNP=ids$SNP,ref_frequency=fr,missingness=miss,ref_pass=good)
    fwrite(ref,file.path(out,"reference_QC.csv"))
    cat("Reference participants: ",nrow(G),"\n")
    fit_one <- function(z,R,n,label) {
      cat("SuSiE: ",label," n=",n," variants=",length(z),"\n",sep="")
      fit <- tryCatch(susieR::susie_rss(z=z,R=R,n=n,L=10,coverage=.95,min_abs_corr=.5,
                 estimate_residual_variance=FALSE,max_iter=1000),error=function(e)e)
      saveRDS(fit,file.path(out,paste0(label,"_fit.rds")))
      if(inherits(fit,"error")) {
        statuses[[label]] <<- data.table(fit=label,n=n,converged=FALSE,nCS=NA_integer_,error=conditionMessage(fit))
        return(NULL)
      }
      ncs <- length(fit$sets$cs)
      statuses[[label]] <<- data.table(fit=label,n=n,converged=isTRUE(fit$converged),nCS=ncs,error="")
      fwrite(data.table(SNP=names(z),PIP=fit$pip),file.path(out,paste0(label,"_PIP.csv")))
      capture.output(str(fit$sets),file=file.path(out,paste0(label,"_credible_sets.txt")))
      if(!isTRUE(fit$converged)) return(NULL)
      fit
    }
    cat("STEP 3/4: Frequency alignment, LD and colocalization\n")
    for(g in names(all_data)) {
      d <- merge(all_data[[g]],ref,by="SNP",all.x=TRUE)
      d[,reason:="keep"]
      d[is.na(ref_pass)|!ref_pass,reason:="reference_MAF_or_missingness"]
      d[reason=="keep" & (abs(xf-ref_frequency)>.15|abs(yf-ref_frequency)>.15),reason:="frequency_mismatch"]
      d[reason=="keep" & palindrome & !(abs(xf-ref_frequency)<=.10 & abs(yf-ref_frequency)<=.10 &
        abs(1-xf-ref_frequency)-abs(xf-ref_frequency)>.10 & abs(1-yf-ref_frequency)-abs(yf-ref_frequency)>.10),reason:="palindrome_frequency_ambiguous"]
      fwrite(d,file.path(out,paste0(g,"_frequency_audit.csv")))
      d <- d[reason=="keep"]; setorder(d,BP,SNP)
      if(nrow(d)<50) stop(g,": fewer than 50 variants after QC")
      if(!all(unname(sentinels) %in% d$SNP)) stop(g,": a key chr15 sentinel failed QC; review audit before fitting")
      stopifnot(length(unique(d$N))==1)
      R <- cor(G[,d$SNP,drop=FALSE]); dimnames(R)<-list(d$SNP,d$SNP)
      stopifnot(all(is.finite(R)),max(abs(R-t(R)))<1e-8,max(abs(diag(R)-1))<1e-8)
      saveRDS(R,file.path(out,paste0(g,"_LD.rds")))
      fwrite(d,file.path(out,paste0(g,"_analysis_input.csv")))
      counts[[g]][, `:=`(final_variants=nrow(d),protein_N=unique(d$N),PTSD_NTOT_median=median(d$NTOT),PTSD_NEFF_median=median(d$NEFF))]
      zx<-setNames(d$xb/d$SE,d$SNP); zy<-setNames(d$yb/d$YSE,d$SNP)
      # Diagnostics are saved for review, not used to silently delete signals.
      diagnostics <- lapply(list(protein=list(z=zx,n=unique(d$N)),PTSD=list(z=zy,n=median(d$NTOT))),function(a)
        tryCatch(susieR::estimate_s_rss(z=a$z,R=R,n=a$n),error=function(e)conditionMessage(e)))
      capture.output(print(diagnostics),file=file.path(out,paste0(g,"_LD_consistency.txt")))
      dx <- list(beta=d$xb,varbeta=d$SE^2,snp=d$SNP,position=d$BP,type="quant",N=unique(d$N),MAF=pmin(d$xf,1-d$xf))
      dy <- list(beta=d$yb,varbeta=d$YSE^2,snp=d$SNP,position=d$BP,type="cc",N=median(d$NTOT),s=median(d$NCASE/d$NTOT))
      for(pr in c(1e-7,5e-7,1e-6,5e-6,1e-5,5e-5,1e-4)) {
        a <- coloc::coloc.abf(dx,dy,p1=1e-4,p2=1e-4,p12=pr)
        abfs[[paste(g,pr)]] <- cbind(data.table(protein=g,p12=pr),as.data.table(as.list(a$summary)))
      }
      fx <- fit_one(zx,R,unique(d$N),paste0(g,"_protein"))
      for(mode in c("NTOT","NEFF")) {
        fy <- fit_one(zy,R,median(d[[mode]]),paste0(g,"_PTSD_",mode))
        if(!is.null(fx)&&!is.null(fy)&&length(fx$sets$cs)>0&&length(fy$sets$cs)>0) {
          for(pr in c(1e-7,5e-7,1e-6,5e-6,1e-5,5e-5,1e-4)) {
            cs <- tryCatch(coloc::coloc.susie(fx,fy,p1=1e-4,p2=1e-4,p12=pr),error=function(e)e)
            saveRDS(cs,file.path(out,paste0(g,"_coloc_",mode,"_",pr,".rds")))
            if(!inherits(cs,"error")&&!is.null(cs$summary)) susies[[paste(g,mode,pr)]] <- cbind(data.table(protein=g,N_mode=mode,p12=pr),as.data.table(cs$summary))
            else cat("Coloc unavailable: ",g," ",mode," ",pr,"\n")
          }
        } else cat("No interpretable SuSiE coloc for ",g,"/",mode,": inspect convergence and credible sets.\n",sep="")
      }
    }
    cat("STEP 4/4: Save review outputs\n")
    fwrite(rbindlist(counts,fill=TRUE),file.path(out,"summary_counts.csv"))
    fwrite(rbindlist(statuses,fill=TRUE),file.path(out,"fit_status.csv"))
    fwrite(rbindlist(abfs,fill=TRUE),file.path(out,"coloc_ABF_prior_grid.csv"))
    if(length(susies)) fwrite(rbindlist(susies,fill=TRUE),file.path(out,"coloc_SuSiE_prior_grid.csv"))
    writeLines("Computation finished; scientific interpretation requires review of QC, convergence, credible sets and LD diagnostics.",file.path(out,"SUCCESS.txt"))
    print(rbindlist(counts,fill=TRUE));print(rbindlist(statuses,fill=TRUE))
    cat("Send summary_counts.csv, fit_status.csv, coloc_*_prior_grid.csv, *_LD_consistency.txt and RUN.log.\n")
    cat("Output directory: ",out,"\n",sep="")
  },error=function(e){writeLines(conditionMessage(e),file.path(out,"FAILED.txt"));message(conditionMessage(e))},
  finally={capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))})
})()
