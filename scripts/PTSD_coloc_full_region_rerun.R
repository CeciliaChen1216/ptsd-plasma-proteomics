# Dense regional colocalization for all ten PTSD-prioritized proteins.
# Revision 2026.10.20.1. Extends the reviewed FES/FURIN implementation.
# Current entry: source('00_START.R'); ptsd_run('coloc')
# Optional before source: options(ptsd.coloc.full_region.resume_dir='...')
# Existing analysis files and released manuscript snapshots are never overwritten.

.pcfr_version <- '2026.10.20.1'
.pcfr_files <- vapply(sys.frames(), function(e) {
  if(is.character(e$ofile) && length(e$ofile)==1L) e$ofile else NA_character_
}, character(1))
.pcfr_files <- .pcfr_files[!is.na(.pcfr_files)]
.pcfr_source_file <- if(length(.pcfr_files)) tail(.pcfr_files,1L) else ''

.pcfr_targets <- function() {
  data.frame(
    protein=c('FES','FURIN','AKT3','CD40','CGREF1','KHK','SIRPA','SNX18','UBE2L6','CD101'),
    chr=c(15L,15L,1L,20L,2L,2L,20L,5L,11L,1L),
    sentinel=c('rs1894401','rs2071410','rs4515770','rs4810485','rs2384572',
               'rs2304681','rs6136377','rs2548611','rs28362950','rs12130298'),
    sentinel_hg19=c(91429042L,91420940L,243686913L,44747947L,27324751L,
                    27315252L,1896288L,53818895L,57371918L,117544331L),
    filename=c('discovery_chr15_FES_P07332_OID21207_v1_Oncology.gz',
      'discovery_chr15_FURIN_P09958_OID21514_v1_Oncology.gz',
      'discovery_chr1_AKT3_Q9Y243_OID21197_v1_Oncology.gz',
      'discovery_chr20_CD40_P25942_OID20724_v1_Inflammation.gz',
      'discovery_chr2_CGREF1_Q99674_OID20152_v1_Cardiometabolic.gz',
      'discovery_chr2_KHK_P50053_OID30241_v1_Cardiometabolic_II.gz',
      'discovery_chr20_SIRPA_P78324_OID20304_v1_Cardiometabolic.gz',
      'discovery_chr5_SNX18_Q96RF0_OID31476_v1_Oncology_II.gz',
      'discovery_chr11_UBE2L6_O14933_OID30321_v1_Cardiometabolic_II.gz',
      'discovery_chr1_CD101_Q93033_OID31480_v1_Oncology_II.gz'),
    stringsAsFactors=FALSE)
}

.pcfr_config <- function() {
  P <- getOption('ptsd.current.paths',list())
  project <- getOption('ptsd.project_root',if(!is.null(P$project_root)) P$project_root else stop('Source 00_START.R or set options(ptsd.project_root=...) first'))
  pick <- function(key,default) if(!is.null(P[[key]]) && nzchar(P[[key]])) P[[key]] else default
  parent <- dirname(project)
  cfg <- list(project_root=project,
    protein_root=pick('protein_raw_base',file.path(parent,'Longcovid/data/raw/pqtl_ukb_ppp/cis_regions')),
    source_casecontrol=pick('source_casecontrol',file.path(project,'data_source_check_20261016/eur_ptsdcasecontrol_pcs_v4_aug3_2021.vcf.gz')),
    pqtl_rds=pick('pqtl_rds',file.path(project,'pQTL/ukbppp_cis_pqtl.rds')),
    plink2=pick('plink2',file.path(parent,'Longcovid/tools',if(.Platform$OS.type=='windows') 'plink2.exe' else 'plink2')),
    ld_prefix=pick('ld_prefix',file.path(parent,'reference/1000G_v3/EUR')),
    output_parent=pick('result_dir',file.path(project,'results')),
    reviewed_snapshot=pick('snapshot_dir',''),
    genes=.pcfr_targets()$protein,chunk_rows=100000L,threads=2L,
    window_bp=500000L,expected_reference_n=503L,info_min=.8,reference_maf_min=.01,
    missingness_max=.05,frequency_difference_max=.15,palindrome_direct_max=.10,
    palindrome_margin=.10,L=10L,coverage=.95,min_abs_corr=.5,max_iter=1000L,
    p1=1e-4,p2=1e-4,p12_reference=1e-5,
    p12_grid=c(1e-7,5e-7,1e-6,5e-6,1e-5,5e-5,1e-4))
  over <- getOption('ptsd.coloc.full_region.overrides',list())
  if(!is.list(over) || (length(over) && (is.null(names(over)) || any(!nzchar(names(over))))))
    stop('ptsd.coloc.full_region.overrides must be a named list')
  unknown <- setdiff(names(over),names(cfg));if(length(unknown)) stop('Unknown overrides: ',paste(unknown,collapse=', '))
  cfg[names(over)] <- over
  if(!length(cfg$genes) || anyDuplicated(cfg$genes) || !all(cfg$genes %in% .pcfr_targets()$protein))
    stop('genes must be a unique subset of the ten prespecified proteins')
  cfg
}

.pcfr_save <- function(x,path) {
  tmp <- paste0(path,'.tmp');saveRDS(x,tmp)
  if(file.exists(path)) unlink(path)
  if(!file.rename(tmp,path)) stop('Cannot finish checkpoint: ',path)
  invisible(path)
}

.pcfr_write <- function(x,path) {
  if(is.null(x) || !ncol(x)) return(invisible(NULL))
  data.table::fwrite(x,path,na='NA')
}

.pcfr_chr <- function(x) suppressWarnings(as.integer(sub('^chr','',as.character(x),ignore.case=TRUE)))
.pcfr_pair <- function(a,b) paste(pmin(toupper(a),toupper(b)),pmax(toupper(a),toupper(b)),sep='/')

.pcfr_fail_if <- function(d,pass,reason) {
  bad <- d$reason=='keep' & (is.na(pass) | !pass)
  data.table::set(d,which(bad),'reason',reason)
  invisible(d)
}

.pcfr_reference <- function(cfg,targets) {
  cat('Reading reference BIM and locating GRCh37 sentinels...\n')
  b <- data.table::fread(paste0(cfg$ld_prefix,'.bim'),header=FALSE,
    col.names=c('CHR','SNP','CM','BP','B1','B2'))
  b[,CHR:=.pcfr_chr(CHR)];b[,`:=`(B1=toupper(B1),B2=toupper(B2))]
  dup <- unique(b$SNP[duplicated(b$SNP)|duplicated(b$SNP,fromLast=TRUE)])
  a <- b[SNP %in% targets$sentinel & !SNP %in% dup]
  if(nrow(a)!=nrow(targets) || anyDuplicated(a$SNP)) stop('A sentinel is absent or duplicated in the reference BIM')
  a <- a[match(targets$sentinel,SNP)]
  if(!all(a$CHR==targets$chr)) stop('Sentinel chromosome differs from the prespecified chromosome')
  if(!all(a$BP==targets$sentinel_hg19)) stop('Reference sentinel positions differ from verified hg19 coordinates; check the reference genome build')
  targets$center_hg19 <- a$BP
  targets$from_bp <- pmax(1L,targets$center_hg19-cfg$window_bp)
  targets$to_bp <- targets$center_hg19+cfg$window_bp
  pieces <- lapply(seq_len(nrow(targets)),function(i) {
    ch <- targets$chr[i];lo <- targets$from_bp[i];hi <- targets$to_bp[i]
    b[CHR==ch & BP>=lo & BP<=hi]
  })
  z <- unique(data.table::rbindlist(pieces))
  z <- z[!SNP %in% dup]
  z[,pair:=.pcfr_pair(B1,B2)]
  list(targets=targets,bim=z,duplicate_ids=dup)
}

.pcfr_read_outcome <- function(path,ids,out,chunk_rows=100000L) {
  # The verified case-control source is used, not the old chromosome-15 subset.
  con <- if(grepl('[.]gz$',path,ignore.case=TRUE)) gzfile(path,'rt') else file(path,'rt')
  on.exit(close(con),add=TRUE)
  header <- NULL;meta <- character()
  repeat {
    line <- readLines(con,n=1L,warn=FALSE)
    if(!length(line)) stop('No case-control header found')
    line <- sub('^\ufeff','',line)
    if(!nzchar(trimws(line))) next
    if(startsWith(line,'##')) {meta<-c(meta,line);next}
    header <- strsplit(sub('^#','',line),'\t',fixed=TRUE)[[1]]
    break
  }
  req <- c('CHROM','ID','POS','A1','A2','FREQ','BETA','SE','P','NCASE','NCON','NEFF','NTOT')
  missing <- setdiff(req,header)
  if(length(missing)) stop('Case-control source lacks columns: ',paste(missing,collapse=', '),
    '. Do not substitute the combined PTSD Z-score file or a .ma file without sample-size fields.')
  writeLines(c(meta,paste(header,collapse='\t')),file.path(out,'PTSD_source_header.txt'),useBytes=TRUE)
  kept <- list();scanned<-0L;nk<-0L
  repeat {
    lines <- readLines(con,n=chunk_rows,warn=FALSE)
    if(!length(lines)) break
    lines <- lines[nzchar(trimws(lines))]
    if(!length(lines)) next
    d <- data.table::fread(text=paste0(paste(lines,collapse='\n'),'\n'),sep='\t',header=FALSE,
      col.names=header,colClasses=list(character=match(c('ID','A1','A2','P'),header)),
      na.strings=c('NA','NaN','.'))
    scanned<-scanned+nrow(d)
    hit <- d[ID %in% ids]
    if(nrow(hit)) {nk<-nk+1L;kept[[nk]]<-hit}
    if(scanned %% (10L*chunk_rows)==0L) cat('PTSD: scanned ',format(scanned,big.mark=','),' rows\n',sep='')
  }
  if(!length(kept)) stop('No target reference SNPs were found in the case-control source')
  y <- data.table::rbindlist(kept,use.names=TRUE)
  data.table::setnames(y,req,c('YCHR','SNP','YBP','YEA','YOA','YFREQ','YBETA','YSE','YP','NCASE','NCON','NEFF','NTOT'))
  y[,`:=`(YCHR=.pcfr_chr(YCHR),YEA=toupper(YEA),YOA=toupper(YOA))]
  dup <- unique(y$SNP[duplicated(y$SNP)|duplicated(y$SNP,fromLast=TRUE)])
  .pcfr_write(y,file.path(out,'PTSD_target_subset_including_duplicates.csv'))
  result <- list(data=y[!SNP %in% dup],duplicate_ids=dup,scanned=scanned,selected=nrow(y),header=header)
  cat('PTSD: ',format(scanned,big.mark=','),' rows read; ',nrow(y),' target rows; ',length(dup),' duplicate IDs excluded.\n',sep='')
  .pcfr_save(result,file.path(out,'PTSD_target_subset_cache.rds'))
  result
}

.pcfr_read_protein <- function(path,target) {
  req <- c('CHROM','GENPOS','ID','ALLELE0','ALLELE1','A1FREQ','INFO','N','TEST','BETA','SE','LOG10P')
  nm <- names(data.table::fread(path,nrows=0L))
  if(length(setdiff(req,nm))) stop('Protein file lacks columns: ',paste(setdiff(req,nm),collapse=', '))
  p <- withCallingHandlers(data.table::fread(path,select=req,
    colClasses=list(character=c('ID','ALLELE0','ALLELE1','TEST'))),
    warning=function(w) stop('Protein file parsing warning; check file completeness: ',conditionMessage(w),call.=FALSE))
  total <- nrow(p)
  parts <- data.table::tstrsplit(p$ID,':',fixed=TRUE,keep=1:4)
  if(length(parts)<4L) stop('Protein variant IDs do not encode chromosome:position:alleles')
  p[,`:=`(CHR=.pcfr_chr(parts[[1]]),BP=suppressWarnings(as.integer(parts[[2]])),
           ID_A0=toupper(parts[[3]]),ID_A1=toupper(parts[[4]]),
           CHROM=.pcfr_chr(CHROM),ALLELE0=toupper(ALLELE0),ALLELE1=toupper(ALLELE1))]
  invalid <- sum(!is.finite(p$CHR)|!is.finite(p$BP))
  ch <- target$chr;lo <- target$from_bp;hi <- target$to_bp
  p <- p[CHR==ch & BP>=lo & BP<=hi]
  in_window <- nrow(p)
  p <- p[TEST=='ADD'];p[,row_id:=seq_len(.N)]
  p[,pair:=.pcfr_pair(ALLELE0,ALLELE1)]
  list(data=p,total=total,invalid_ids=invalid,in_window=in_window,additive=nrow(p))
}

.pcfr_prepare_region <- function(p,b,y,target,outcome_duplicates=character()) {
  b <- data.table::copy(b)
  ch <- target$chr;lo <- target$from_bp;hi <- target$to_bp
  b <- b[CHR==ch & BP>=lo & BP<=hi]
  # Keep unmatched regional records in the audit, rather than only an inner join.
  d <- merge(p,b[,.(CHR,BP,pair,SNP,CM,B1,B2)],by=c('CHR','BP','pair'),all.x=TRUE,allow.cartesian=TRUE)
  d <- merge(d,y,by='SNP',all.x=TRUE,allow.cartesian=TRUE)
  d[,reason:='keep']
  .pcfr_fail_if(d,is.finite(d$BP)&d$CHROM==d$CHR&d$pair==.pcfr_pair(d$ID_A0,d$ID_A1),'protein_ID_header_mismatch')
  .pcfr_fail_if(d,d$ALLELE0 %in% c('A','C','G','T') & d$ALLELE1 %in% c('A','C','G','T') & d$ALLELE0!=d$ALLELE1,'not_biallelic_SNV')
  .pcfr_fail_if(d,!is.na(d$SNP)&nzchar(d$SNP),'no_reference_position_allele_match')
  duplicate <- duplicated(d$row_id)|duplicated(d$row_id,fromLast=TRUE)|
    (!is.na(d$SNP)&(duplicated(d$SNP)|duplicated(d$SNP,fromLast=TRUE)))
  .pcfr_fail_if(d,!duplicate,'duplicate_mapping')
  .pcfr_fail_if(d,!d$SNP %in% outcome_duplicates,'outcome_duplicate_ID')
  .pcfr_fail_if(d,is.finite(d$YBP),'outcome_SNP_absent')
  .pcfr_fail_if(d,d$CHR==d$YCHR & d$BP==d$YBP &
    ((d$YEA==d$B1 & d$YOA==d$B2)|(d$YEA==d$B2 & d$YOA==d$B1)),'outcome_pair_or_position_mismatch')
  valid <- is.finite(d$BETA)&is.finite(d$SE)&d$SE>0&is.finite(d$YBETA)&is.finite(d$YSE)&d$YSE>0&
    is.finite(d$INFO)&is.finite(d$N)&d$N>0&is.finite(d$NTOT)&d$NTOT>0&is.finite(d$NEFF)&d$NEFF>0&
    is.finite(d$NCASE)&d$NCASE>0&is.finite(d$NCON)&d$NCON>0&
    is.finite(d$A1FREQ)&d$A1FREQ>0&d$A1FREQ<1&is.finite(d$YFREQ)&d$YFREQ>0&d$YFREQ<1
  .pcfr_fail_if(d,valid,'invalid_effect_SE_frequency_or_N')
  .pcfr_fail_if(d,abs(d$NTOT-d$NCASE-d$NCON)<=1 & d$NEFF<=d$NTOT,'outcome_sample_size_inconsistent')
  # INFO threshold is set by the caller so it can be recorded in the run config.
  d[,`:=`(xb=ifelse(ALLELE1==B1,BETA,-BETA),yb=ifelse(YEA==B1,YBETA,-YBETA),
    xf=ifelse(ALLELE1==B1,A1FREQ,1-A1FREQ),yf=ifelse(YEA==B1,YFREQ,1-YFREQ),
    palindrome=pair %in% c('A/T','C/G'))]
  d
}

.pcfr_genotypes <- function(d,cfg,out) {
  ids <- d[,.(SNP,B1)]
  if(anyDuplicated(ids$SNP)) stop('Duplicate genotype export IDs')
  .pcfr_write(ids[,.(SNP)],file.path(out,'extract_with_header.csv'))
  data.table::fwrite(ids[,.(SNP)],file.path(out,'extract.txt'),col.names=FALSE)
  data.table::fwrite(ids,file.path(out,'counted_alleles.txt'),sep='\t',col.names=FALSE)
  prefix <- file.path(out,'reference')
  args <- c('--bfile',shQuote(cfg$ld_prefix),'--extract',shQuote(file.path(out,'extract.txt')),
    '--export','A','--export-allele',shQuote(file.path(out,'counted_alleles.txt')),
    '--threads',as.character(cfg$threads),'--out',shQuote(prefix))
  writeLines(paste(shQuote(cfg$plink2),paste(args,collapse=' ')),file.path(out,'PLINK.command.txt'))
  code <- system2(cfg$plink2,args,stdout=file.path(out,'PLINK_console.txt'),stderr=file.path(out,'PLINK_stderr.txt'))
  if(code!=0 || !file.exists(paste0(prefix,'.raw'))) stop('PLINK export failed; inspect PLINK logs')
  raw <- data.table::fread(paste0(prefix,'.raw'),na.strings=c('NA','nan'))
  ans <- .pcfr_decode_genotypes(raw,ids,cfg)
  .pcfr_write(ans$ref,file.path(out,'reference_QC.csv'))
  ans
}

.pcfr_decode_genotypes <- function(raw,ids,cfg) {
  columns <- paste(ids$SNP,ids$B1,sep='_')
  if(!all(columns %in% names(raw))) stop('PLINK counted-allele column names do not match the requested B1 alleles')
  if(nrow(raw)!=cfg$expected_reference_n) stop('Exported reference N differs from expected reference N')
  G <- as.matrix(raw[,..columns]);storage.mode(G)<-'double';colnames(G)<-ids$SNP
  if(any(G<0|G>2,na.rm=TRUE)) stop('Unexpected genotype dosage outside 0 to 2')
  miss <- colMeans(is.na(G));fr <- colMeans(G,na.rm=TRUE)/2
  for(j in seq_len(ncol(G))) G[is.na(G[,j]),j]<-2*fr[j]
  good <- is.finite(fr)&pmin(fr,1-fr)>=cfg$reference_maf_min&miss<=cfg$missingness_max
  ref <- data.table::data.table(SNP=ids$SNP,ref_frequency=unname(fr),missingness=unname(miss),ref_pass=unname(good))
  list(G=G[,good,drop=FALSE],ref=ref)
}

.pcfr_frequency <- function(d,ref,cfg) {
  z <- merge(d,ref,by='SNP',all.x=TRUE)
  z[,reason:='keep']
  .pcfr_fail_if(z,z$ref_pass,'reference_MAF_or_missingness')
  .pcfr_fail_if(z,abs(z$xf-z$ref_frequency)<=cfg$frequency_difference_max &
    abs(z$yf-z$ref_frequency)<=cfg$frequency_difference_max,'frequency_mismatch')
  pass <- !z$palindrome | (abs(z$xf-z$ref_frequency)<=cfg$palindrome_direct_max &
    abs(z$yf-z$ref_frequency)<=cfg$palindrome_direct_max &
    abs(1-z$xf-z$ref_frequency)-abs(z$xf-z$ref_frequency)>cfg$palindrome_margin &
    abs(1-z$yf-z$ref_frequency)-abs(z$yf-z$ref_frequency)>cfg$palindrome_margin)
  .pcfr_fail_if(z,pass,'palindrome_frequency_ambiguous')
  z
}

.pcfr_sentinel <- function(d,target,catalog,out) {
  q <- catalog[protein==target$protein & SNP==target$sentinel]
  a <- d[SNP==target$sentinel]
  if(nrow(q)!=1L || nrow(a)!=1L) stop('The direct sentinel is absent, duplicated or failed QC')
  if(.pcfr_pair(q$effect_allele,q$other_allele)!=a$pair) stop('Sentinel allele pair differs from the original instrument catalog')
  observed <- if(toupper(q$effect_allele)==a$B1) a$xb else -a$xb
  match_beta <- abs(observed-q$beta)<=1e-8+1e-5*abs(q$beta)
  match_se <- abs(a$SE-q$se)<=1e-8+1e-5*abs(q$se)
  audit <- data.table::data.table(protein=target$protein,SNP=target$sentinel,
    source_effect_allele=q$effect_allele,source_other_allele=q$other_allele,
    catalog_beta=q$beta,raw_beta_aligned_to_catalog=observed,catalog_se=q$se,raw_se=a$SE,
    beta_matches=match_beta,se_matches=match_se)
  .pcfr_write(audit,file.path(out,'sentinel_catalog_check.csv'))
  if(!match_beta || !match_se) stop('Protein sentinel beta/SE do not reproduce the original instrument catalog; check assay/file identity')
  if(target$protein %in% c('FES','FURIN') && !all(c('rs1894401','rs2071410') %in% d$SNP))
    stop('A required chromosome-15 sentinel failed QC')
  invisible(audit)
}

.pcfr_datasets <- function(d) {
  np <- unique(d$N)
  if(length(np)!=1L || !is.finite(np) || np<=0) stop('Variable or invalid protein N: review rather than silently use a scalar median')
  list(x=list(beta=d$xb,varbeta=d$SE^2,snp=d$SNP,position=d$BP,type='quant',N=np,MAF=pmin(d$xf,1-d$xf)),
    y=list(beta=d$yb,varbeta=d$YSE^2,snp=d$SNP,position=d$BP,type='cc',N=median(d$NTOT),s=median(d$NCASE/d$NTOT)),
    zx=stats::setNames(d$xb/d$SE,d$SNP),zy=stats::setNames(d$yb/d$YSE,d$SNP))
}

.pcfr_fit <- function(z,R,n,label,cfg,out) {
  cat('SuSiE ',label,': ',length(z),' SNPs; N=',n,'\n',sep='')
  warnings <- character()
  fit <- withCallingHandlers(tryCatch(susieR::susie_rss(z=z,R=R,n=n,L=cfg$L,
    coverage=cfg$coverage,min_abs_corr=cfg$min_abs_corr,estimate_residual_variance=FALSE,
    max_iter=cfg$max_iter),error=function(e)e),warning=function(w) warnings<<-c(warnings,conditionMessage(w)))
  .pcfr_save(fit,file.path(out,paste0(label,'_fit.rds')))
  error <- if(inherits(fit,'error')) conditionMessage(fit) else ''
  converged <- !nzchar(error) && isTRUE(fit$converged)
  ncs <- if(nzchar(error)) NA_integer_ else length(fit$sets$cs)
  status <- data.table::data.table(fit=label,n=n,converged=converged,nCS=ncs,error=error,
    warnings=paste(unique(warnings),collapse='; '))
  .pcfr_write(status,file.path(out,paste0(label,'_status.csv')))
  if(!nzchar(error)) {
    .pcfr_write(data.table::data.table(SNP=names(z),PIP=fit$pip),file.path(out,paste0(label,'_PIP.csv')))
    capture.output(str(fit$sets),file=file.path(out,paste0(label,'_credible_sets.txt')))
    cs_names <- names(fit$sets$cs)
    if(is.null(cs_names)) cs_names <- paste0('CS',seq_along(fit$sets$cs))
    members <- lapply(seq_along(fit$sets$cs),function(j) data.table::data.table(
      credible_set=cs_names[j],SNP=names(z)[fit$sets$cs[[j]]],index=fit$sets$cs[[j]]))
    if(length(members)) .pcfr_write(data.table::rbindlist(members),file.path(out,paste0(label,'_credible_set_members.csv')))
  }
  list(fit=if(converged) fit else NULL,status=status)
}

.pcfr_control <- function(target,d,susie,abf,cfg) {
  if(!target$protein %in% c('FES','FURIN')) return(NULL)
  expected <- c(FES=2516L,FURIN=2517L)[target$protein]
  z <- data.table::data.table(protein=target$protein,metric='final_variant_count',N_mode='both',
    observed=nrow(d),reviewed=unname(expected),difference=nrow(d)-expected,matches=nrow(d)==expected)
  if(nzchar(cfg$reviewed_snapshot)) {
    for(method in c('SuSiE','ABF')) {
      f <- file.path(cfg$reviewed_snapshot,paste0('coloc_',method,'_prior_grid.csv'))
      if(!file.exists(f)) next
      old <- data.table::fread(f);new <- if(method=='SuSiE') susie else abf
      if(is.null(new) || !nrow(new)) next
      modes <- if(method=='SuSiE') c('NTOT','NEFF') else 'NTOT'
      for(mode in modes) {
        a <- old[protein==target$protein & abs(p12-cfg$p12_reference)<1e-15]
        if(method=='SuSiE') a<-a[N_mode==mode]
        v <- new[abs(p12-cfg$p12_reference)<1e-15]
        if(method=='SuSiE') v<-v[N_mode==mode]
        if(!nrow(a)||!nrow(v)||!any(is.finite(a$PP.H4.abf))||!any(is.finite(v$PP.H4.abf))) next
        before<-max(a$PP.H4.abf,na.rm=TRUE);after<-max(v$PP.H4.abf,na.rm=TRUE)
        z<-data.table::rbindlist(list(z,data.table::data.table(protein=target$protein,
          metric=paste0(method,'_maximum_pair_PP_H4_reference_prior'),N_mode=mode,
          observed=after,reviewed=before,difference=after-before,matches=abs(after-before)<=1e-6)),fill=TRUE)
      }
    }
  }
  z
}

.pcfr_region <- function(target,b,y,ydup,catalog,cfg,out) {
  p <- .pcfr_read_protein(target$path,target)
  if(!nrow(p$data)) stop('No additive protein records in the hg19 window parsed from variant IDs')
  d <- .pcfr_prepare_region(p$data,b,y,target,ydup)
  .pcfr_fail_if(d,d$INFO>=cfg$info_min,'protein_INFO_below_threshold')
  .pcfr_write(d,file.path(out,'preLD_audit.csv'))
  .pcfr_write(d[,.N,by=reason],file.path(out,'preLD_reason_counts.csv'))
  counts <- data.table::data.table(protein=target$protein,sentinel=target$sentinel,chr=target$chr,
    center_hg19=target$center_hg19,from_bp=target$from_bp,to_bp=target$to_bp,
    protein_chromosome_rows=p$total,unparseable_ID_rows=p$invalid_ids,
    protein_window_rows=p$in_window,additive_region_rows=p$additive,
    reference_position_pair_matches=sum(!is.na(d$SNP)),outcome_present=sum(is.finite(d$YBP)),
    preLD_pass=sum(d$reason=='keep'))
  .pcfr_write(counts,file.path(out,'summary_counts.csv'))
  d <- d[reason=='keep']
  if(nrow(d)<50L) stop('Fewer than 50 aligned variants before LD QC; inspect the full-region audit')
  gt <- .pcfr_genotypes(d,cfg,out)
  q <- .pcfr_frequency(d,gt$ref,cfg)
  .pcfr_write(q,file.path(out,'frequency_audit.csv'))
  .pcfr_write(q[,.N,by=reason],file.path(out,'frequency_reason_counts.csv'))
  d <- q[reason=='keep'];data.table::setorder(d,BP,SNP)
  counts[,`:=`(final_variants=nrow(d),reference_participants=nrow(gt$G))]
  .pcfr_write(counts,file.path(out,'summary_counts.csv'))
  if(nrow(d)<50L) stop('Fewer than 50 variants after LD/frequency QC; not an interpretable dense-region fit')
  .pcfr_sentinel(d,target,catalog,out)
  ds <- .pcfr_datasets(d)
  R <- stats::cor(gt$G[,d$SNP,drop=FALSE]);dimnames(R)<-list(d$SNP,d$SNP)
  stopifnot(all(is.finite(R)),max(abs(R-t(R)))<1e-8,max(abs(diag(R)-1))<1e-8,
    identical(colnames(R),ds$x$snp),identical(names(ds$zx),ds$x$snp))
  .pcfr_save(R,file.path(out,'LD.rds'))
  .pcfr_write(d,file.path(out,'analysis_input.csv'))
  counts[,`:=`(protein_N=ds$x$N,PTSD_NTOT_median=median(d$NTOT),PTSD_NEFF_median=median(d$NEFF),
    max_protein_LOG10P=max(d$LOG10P,na.rm=TRUE),min_PTSD_P_from_Z=2*pnorm(-max(abs(ds$zy))))]
  .pcfr_write(counts,file.path(out,'summary_counts.csv'))
  diagnostic <- lapply(list(protein=list(z=ds$zx,n=ds$x$N),PTSD=list(z=ds$zy,n=median(d$NTOT))),function(a)
    tryCatch(susieR::estimate_s_rss(z=a$z,R=R,n=a$n),error=function(e)conditionMessage(e)))
  capture.output(print(diagnostic),file=file.path(out,paste0(target$protein,'_LD_consistency.txt')))
  fits <- list();abfs<-list();pairs<-list();pair_status<-list()
  # Full-region ABF is reported even if SuSiE does not identify credible sets.
  for(pr in cfg$p12_grid) {
    a <- tryCatch(coloc::coloc.abf(ds$x,ds$y,p1=cfg$p1,p2=cfg$p2,p12=pr),error=function(e)e)
    if(inherits(a,'error')) arow<-data.table::data.table(protein=target$protein,p12=pr,N_mode='NTOT',
      error=conditionMessage(a),PP.H0.abf=NA_real_,PP.H1.abf=NA_real_,PP.H2.abf=NA_real_,PP.H3.abf=NA_real_,PP.H4.abf=NA_real_)
    else arow<-cbind(data.table::data.table(protein=target$protein,p12=pr,N_mode='NTOT',error=''),data.table::as.data.table(as.list(a$summary)))
    abfs[[as.character(pr)]]<-arow
    .pcfr_write(data.table::rbindlist(abfs,fill=TRUE),file.path(out,'coloc_ABF_prior_grid.csv'))
  }
  # This execution order reproduces the reviewed CD101 completion. It changes
  # neither the SNP set nor the model and records an unrun protein fit as NA.
  outcome_cache <- list()
  if(target$protein=='CD101') {
    timed_fit <- function(z,n,label) R.utils::withTimeout(
      .pcfr_fit(z,R,n,label,cfg,out),timeout=900,onTimeout='error')
    for(mode in c('NTOT','NEFF'))
      outcome_cache[[mode]] <- timed_fit(ds$zy,median(d[[mode]]),paste0(target$protein,'_PTSD_',mode))
    no_outcome_cs <- all(vapply(outcome_cache,function(x)
      !is.null(x$fit) && length(x$fit$sets$cs)==0L,logical(1)))
    if(no_outcome_cs) {
      fx <- list(fit=NULL,status=data.table::data.table(fit=paste0(target$protein,'_protein'),
        n=ds$x$N,converged=NA,nCS=NA_integer_,error='',warnings='',
        execution='not_run_no_PTSD_credible_set'))
      .pcfr_write(fx$status,file.path(out,paste0(target$protein,'_protein_status.csv')))
    } else fx <- timed_fit(ds$zx,ds$x$N,paste0(target$protein,'_protein'))
  } else fx <- .pcfr_fit(ds$zx,R,ds$x$N,paste0(target$protein,'_protein'),cfg,out)
  fits[[1]] <- fx$status
  for(mode in c('NTOT','NEFF')) {
    fy <- if(length(outcome_cache)) outcome_cache[[mode]] else
      .pcfr_fit(ds$zy,R,median(d[[mode]]),paste0(target$protein,'_PTSD_',mode),cfg,out)
    fits[[length(fits)+1L]]<-fy$status
    .pcfr_write(data.table::rbindlist(fits,fill=TRUE),file.path(out,'fit_status.csv'))
    why <- if(target$protein=='CD101') .pcfr_completion_pair_status(fx,fy) else
      if(is.null(fx$fit)||is.null(fy$fit)) 'fit_error_or_nonconvergence' else
      if(!length(fx$fit$sets$cs) && !length(fy$fit$sets$cs)) 'no_protein_or_PTSD_credible_set' else
      if(!length(fx$fit$sets$cs)) 'no_protein_credible_set' else
      if(!length(fy$fit$sets$cs)) 'no_PTSD_credible_set' else ''
    for(pr in cfg$p12_grid) {
      label<-paste(mode,pr,sep='_');err<-'';npairs<-0L
      if(nzchar(why)) status<-why else {
        cs<-tryCatch(coloc::coloc.susie(fx$fit,fy$fit,p1=cfg$p1,p2=cfg$p2,p12=pr),error=function(e)e)
        .pcfr_save(cs,file.path(out,paste0('coloc_',label,'.rds')))
        if(inherits(cs,'error')) {status<-'coloc_error';err<-conditionMessage(cs)} else
          if(is.null(cs$summary)||!nrow(cs$summary)) status<-'no_signal_pairs_returned' else {
            status<-'estimated';npairs<-nrow(cs$summary)
            pairs[[label]]<-cbind(data.table::data.table(protein=target$protein,N_mode=mode,p12=pr),data.table::as.data.table(cs$summary))
          }
      }
      pair_status[[label]]<-data.table::data.table(protein=target$protein,N_mode=mode,p12=pr,status=status,n_pairs=npairs,error=err)
    }
    if(nzchar(why)) cat(target$protein,'/',mode,': ',why,'; ABF remains reported.\n',sep='')
  }
  su <- if(length(pairs)) data.table::rbindlist(pairs,fill=TRUE) else NULL
  ab <- data.table::rbindlist(abfs,fill=TRUE);fs<-data.table::rbindlist(fits,fill=TRUE)
  ps <- data.table::rbindlist(pair_status,fill=TRUE)
  .pcfr_write(su,file.path(out,'coloc_SuSiE_prior_grid.csv'))
  .pcfr_write(ps,file.path(out,'coloc_status.csv'))
  control <- .pcfr_control(target,d,su,ab,cfg)
  .pcfr_write(control,file.path(out,'reviewed_chr15_control_comparison.csv'))
  overview <- data.table::rbindlist(lapply(c('NTOT','NEFF'),function(mode) {
    z<-if(is.null(su)) NULL else su[N_mode==mode & abs(p12-cfg$p12_reference)<1e-15]
    v<-if(is.null(z)||!nrow(z)) NULL else z[which.max(PP.H4.abf)]
    s<-ps[N_mode==mode & abs(p12-cfg$p12_reference)<1e-15]
    data.table::data.table(protein=target$protein,N_mode=mode,n_snps=nrow(d),p12=cfg$p12_reference,
      nCS_protein=fs$nCS[1],nCS_PTSD=fs$nCS[match(paste0(target$protein,'_PTSD_',mode),fs$fit)],
      status=s$status[1],n_pairs=s$n_pairs[1],
      maximum_pair_PP_H4=if(is.null(v)) NA_real_ else v$PP.H4.abf,
      same_pair_PP_H3=if(is.null(v)) NA_real_ else v$PP.H3.abf,
      hit1=if(is.null(v)) NA_character_ else as.character(v$hit1),
      hit2=if(is.null(v)) NA_character_ else as.character(v$hit2))
  }),fill=TRUE)
  result<-list(counts=counts,fit_status=fs,abf=ab,susie=su,pair_status=ps,overview=overview,control=control)
  # Written only after all outputs for this protein have finished.
  .pcfr_save(result,file.path(out,'REGION_COMPLETE.rds'))
  result
}

.pcfr_completion_pair_status <- function(protein,outcome) {
  if(is.null(outcome$fit)) 'fit_error_or_nonconvergence' else
    if(!length(outcome$fit$sets$cs)) 'no_PTSD_credible_set' else
    if(is.null(protein$fit)) 'fit_error_or_nonconvergence' else
    if(!length(protein$fit$sets$cs)) 'no_protein_credible_set' else ''
}

.pcfr_run_status <- function(ans) {
  skipped <- 'execution' %in% names(ans$fit_status) &&
    any(ans$fit_status$execution=='not_run_no_PTSD_credible_set',na.rm=TRUE)
  ab_ok <- !anyNA(ans$abf$error) && all(!nzchar(ans$abf$error))
  pair_ok <- !any(ans$pair_status$status %in% c('coloc_error','no_signal_pairs_returned'))
  if(skipped && all(ans$pair_status$status=='no_PTSD_credible_set') &&
     all(ans$fit_status$converged[!is.na(ans$fit_status$converged)]) && ab_ok && pair_ok)
    'assessment_complete_no_PTSD_credible_set'
  else if(isTRUE(all(ans$fit_status$converged)) && ab_ok && pair_ok) 'completed'
  else 'completed_with_fit_or_coloc_error'
}

.pcfr_aggregate <- function(results,out) {
  for(key in c('counts','fit_status','abf','susie','pair_status','overview','control')) {
    pieces<-lapply(results,function(x)x[[key]]);pieces<-pieces[vapply(pieces,function(x)!is.null(x)&&nrow(x)>0L,logical(1))]
    name<-c(counts='summary_counts.csv',fit_status='fit_status.csv',abf='coloc_ABF_prior_grid.csv',
      susie='coloc_SuSiE_prior_grid.csv',pair_status='coloc_status.csv',overview='reference_prior_overview.csv',
      control='reviewed_chr15_control_comparison.csv')[[key]]
    if(length(pieces)) .pcfr_write(data.table::rbindlist(pieces,fill=TRUE),file.path(out,name))
  }
}

.pcfr_sendback <- function(out) {
  send<-file.path(out,'SEND_BACK');dir.create(send,showWarnings=FALSE)
  small <- c('RUN.log','FAILED.txt','FINISHED.txt','run_status.csv','parameters.txt','input_manifest.csv',
    'target_manifest.csv','package_versions.csv','summary_counts.csv','fit_status.csv','coloc_ABF_prior_grid.csv',
    'coloc_SuSiE_prior_grid.csv','coloc_status.csv','reference_prior_overview.csv',
    'reviewed_chr15_control_comparison.csv','PTSD_source_header.txt','sessionInfo.txt','script_snapshot.R.txt','config_used.txt')
  small<-small[file.exists(file.path(out,small))]
  if(length(small)) file.copy(file.path(out,small),file.path(send,small),overwrite=TRUE)
  regions<-list.dirs(file.path(out,'regions'),recursive=FALSE,full.names=TRUE)
  for(folder in regions) {
    names<-list.files(folder,pattern='LD_consistency[.]txt$|reason_counts[.]csv$|sentinel_catalog_check[.]csv$|FAILED[.]txt$|PLINK.*[.]txt$|[.]log$|summary_counts[.]csv$|_status[.]csv$|_PIP[.]csv$|_credible_set_members[.]csv$|_credible_sets[.]txt$',full.names=TRUE)
    if(length(names)) file.copy(names,file.path(send,paste0(basename(folder),'__',basename(names))),overwrite=TRUE)
  }
  writeLines(c('Review-only bundle; full aligned data, genotypes, LD and fit objects remain in the original result directory.',
    'No-SuSiE-pair statuses distinguish missing credible sets from fitting errors; ABF H0-H4 are retained.',
    'Maximum-pair PP.H4 is descriptive, not a mean or automatic effector assignment.'),file.path(send,'README.txt'))
  old<-getwd();on.exit(setwd(old),add=TRUE);setwd(out)
  if(requireNamespace('zip',quietly=TRUE)) {
    archive<-file.path(out,'SEND_BACK.zip');zip::zipr(archive,files='SEND_BACK')
  } else {
    archive<-file.path(out,'SEND_BACK.tar.gz');utils::tar(archive,files='SEND_BACK',compression='gzip',tar='internal')
  }
  cat('Upload this review bundle: ',normalizePath(archive,winslash='/',mustWork=TRUE),'\n',sep='')
  invisible(archive)
}

ptsd_coloc_full_region <- function() {
  cfg<-.pcfr_config()
  required<-c('data.table','R.utils','coloc','susieR')
  absent<-required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
  if(length(absent)) stop('Missing packages: ',paste(absent,collapse=', '),'. Install them before running; nothing is installed automatically.')
  suppressPackageStartupMessages(library(data.table))
  targets<-data.table::as.data.table(.pcfr_targets());targets<-targets[match(cfg$genes,protein)]
  targets[,path:=file.path(cfg$protein_root,protein,filename)]
  versions<-data.table::data.table(package=c('R',required),version=c(as.character(getRversion()),vapply(required,function(p)as.character(utils::packageVersion(p)),character(1))))
  inputs<-c(PLINK2=cfg$plink2,BED=paste0(cfg$ld_prefix,'.bed'),BIM=paste0(cfg$ld_prefix,'.bim'),
    FAM=paste0(cfg$ld_prefix,'.fam'),PTSD=cfg$source_casecontrol,pQTL_catalog=cfg$pqtl_rds,
    stats::setNames(targets$path,targets$protein))
  resume<-getOption('ptsd.coloc.full_region.resume_dir','')
  if(length(resume)!=1L || is.na(resume)) stop('resume_dir must be a single path or an empty string')
  dir.create(cfg$output_parent,recursive=TRUE,showWarnings=FALSE)
  if(nzchar(resume)) {
    if(!dir.exists(resume)||!file.exists(file.path(resume,'RUN_STATE.rds'))) stop('Explicit resume directory has no RUN_STATE.rds')
    out<-normalizePath(resume,winslash='/',mustWork=TRUE)
  } else {
    out<-tempfile(paste0('coloc_full_region_',format(Sys.time(),'%Y%m%d_%H%M%S'),'_'),tmpdir=cfg$output_parent)
    if(!dir.create(out)) stop('Cannot create output directory')
    out<-normalizePath(out,winslash='/',mustWork=TRUE)
  }
  cat('Output directory: ',out,'\n',sep='')
  log<-file(file.path(out,'RUN.log'),if(nzchar(resume)) 'at' else 'wt')
  sink(log,split=TRUE);sink(log,type='message')
  fatal<-NULL;results<-list();run_status<-list()
  on.exit({
    capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
    sink(type='message');sink();close(log)
    tryCatch(.pcfr_sendback(out),error=function(e)message('Review bundle creation failed; upload SEND_BACK folder manually: ',conditionMessage(e)))
  },add=TRUE)
  tryCatch({
    cat('STEP 1/4: Verify inputs and preserve the original pipeline settings\n')
    if(any(!file.exists(inputs))) stop('Missing inputs: ',paste(inputs[!file.exists(inputs)],collapse='; '))
    info<-file.info(inputs)
    manifest<-data.table::data.table(item=names(inputs),path=unname(inputs),bytes=info$size,mtime=as.character(info$mtime),md5=NA_character_)
    hash<-manifest$item!='BED'
    cat('Recording input checksums (reference BED: size and modification time only)...\n')
    manifest$md5[hash]<-unname(tools::md5sum(manifest$path[hash]))
    if(any(is.na(manifest$md5[hash]))) stop('An input checksum could not be recorded')
    script_md5<-if(nzchar(.pcfr_source_file)&&file.exists(.pcfr_source_file)) unname(tools::md5sum(.pcfr_source_file)) else NA_character_
    signature<-list(version=.pcfr_version,config=cfg,input_manifest=manifest,package_versions=versions,script_md5=script_md5)
    if(nzchar(resume)) {
      saved<-readRDS(file.path(out,'RUN_STATE.rds'))
      if(!identical(saved,signature)) stop('Resume rejected: code, settings, package versions or input fingerprints changed. Clear the resume option and start a new run.')
      cat('Explicit resume accepted; only completed region checkpoints will be reused.\n')
      for(previous in c('FAILED.txt','FINISHED.txt')) if(file.exists(file.path(out,previous)))
        file.rename(file.path(out,previous),file.path(out,paste0('PREVIOUS_',format(Sys.time(),'%Y%m%d_%H%M%S'),'_',previous)))
    } else .pcfr_save(signature,file.path(out,'RUN_STATE.rds'))
    .pcfr_write(manifest,file.path(out,'input_manifest.csv'));.pcfr_write(versions,file.path(out,'package_versions.csv'))
    if(nzchar(.pcfr_source_file)&&file.exists(.pcfr_source_file)) file.copy(.pcfr_source_file,file.path(out,'script_snapshot.R.txt'),overwrite=TRUE)
    writeLines(c(paste('Code revision:',.pcfr_version),
      'Default settings described below; actual settings, including any overrides, are recorded in config_used.txt.',
      'All ten candidates by default; FES/FURIN serve as reviewed dense-region controls.',
      'CD101 outcome fits run first with a 900-second fit time limit; if both converge without a 95% credible set the protein fit is explicitly not run.',
      'Each sentinel +/-500kb in GRCh37; parse protein variant IDs, never use GENPOS as hg19.',
      'No association-P filtering, clumping, SNP-count cap, diagnostic-driven deletion or automatic ridge.',
      'Biallelic A/C/G/T; protein INFO>=0.8; reference MAF>=0.01; genotype missingness<=0.05.',
      'Exact position and allele-pair matching; no automatic strand-complement rescue; duplicated reference/outcome IDs excluded.',
      'Effects and genotypes aligned to BIM B1, explicitly counted with PLINK --export-allele.',
      'Aligned exposure/outcome frequency distance<=0.15; palindromes also require direct distance<=0.10 and alternative-minus-direct distance>0.10.',
      'LD is Pearson correlation after per-variant mean imputation, using the audited 503-participant EUR panel.',
      'SuSiE uses z statistics, L=10, coverage=.95, min_abs_corr=.5, max_iter=1000 and fixed residual variance.',
      'Protein N must be constant in the retained regional file. PTSD primary N=median retained NTOT; sensitivity N=median retained source NEFF.',
      'ABF uses beta/SE, protein MAF and actual N with sdY estimated by coloc; outcome type=cc and source log-odds beta/SE.',
      'ABF is retained when credible sets are absent. SuSiE no-CS, fitting errors, nonconvergence and missing inputs have separate statuses.',
      'p1=p2=1e-4; p12=1e-7,5e-7,1e-6,5e-6,1e-5,5e-5,1e-4; reference p12=1e-5.',
      'All signal pairs and raw credible-set memberships are saved. Overview uses maximum-pair PP.H4, not a mean or automatic causal verdict.',
      'Full original GWAS regional coverage is not claimed: matching and reference QC restrict the available variants.',
      'Fixed minimum of 50 retained SNPs is a technical safeguard inherited from the reviewed implementation, not proof of adequate regional coverage.',
      'Every completed protein is checkpointed. Resuming requires an explicit path and identical code/settings/input fingerprints.'),file.path(out,'parameters.txt'))
    # Record actual overrides as well as the readable description above.
    dput(cfg,file=file.path(out,'config_used.txt'))
    print(system2(cfg$plink2,'--version',stdout=TRUE,stderr=TRUE))
    nfam<-nrow(data.table::fread(paste0(cfg$ld_prefix,'.fam'),header=FALSE))
    if(nfam!=cfg$expected_reference_n) stop('Reference FAM count is not the expected 503 participants')
    ref<-.pcfr_reference(cfg,targets);targets<-ref$targets;b<-ref$bim
    .pcfr_write(targets,file.path(out,'target_manifest.csv'))
    .pcfr_write(b,file.path(out,'regional_reference_BIM.csv'))
    .pcfr_write(data.table::data.table(SNP=ref$duplicate_ids),file.path(out,'excluded_reference_duplicate_IDs.csv'))
    catalog<-data.table::as.data.table(readRDS(cfg$pqtl_rds))
    need<-c('protein','SNP','effect_allele','other_allele','beta','se')
    if(length(setdiff(need,names(catalog)))) stop('Instrument catalog lacks required columns')
    cat('STEP 2/4: Read all target PTSD regions in one pass through the verified case-control source\n')
    cache<-file.path(out,'PTSD_target_subset_cache.rds')
    yy<-if(nzchar(resume)&&file.exists(cache)) readRDS(cache) else .pcfr_read_outcome(cfg$source_casecontrol,unique(b$SNP),out,cfg$chunk_rows)
    dir.create(file.path(out,'regions'),showWarnings=FALSE)
    cat('STEP 3/4: Dense regional QC, reference export, fine-mapping and colocalization\n')
    for(i in seq_len(nrow(targets))) {
      target<-targets[i];g<-target$protein;folder<-file.path(out,'regions',g);dir.create(folder,showWarnings=FALSE)
      cat('\n[',i,'/',nrow(targets),'] ',g,' chr',target$chr,':',target$from_bp,'-',target$to_bp,'\n',sep='')
      checkpoint<-file.path(folder,'REGION_COMPLETE.rds')
      ans<-tryCatch({
        if(nzchar(resume)&&file.exists(checkpoint)) {cat('Using this protein\'s completed checkpoint.\n');readRDS(checkpoint)}
        else .pcfr_region(target,b,yy$data,yy$duplicate_ids,catalog,cfg,folder)
      },error=function(e) {
        writeLines(conditionMessage(e),file.path(folder,'FAILED.txt'))
        cat('REGION NEEDS REVIEW: ',g,': ',conditionMessage(e),'\n',sep='');NULL
      })
      if(is.null(ans)) {
        run_status[[g]]<-data.table::data.table(protein=g,status='region_failed',detail=readLines(file.path(folder,'FAILED.txt'),warn=FALSE)[1])
      } else {
        if(file.exists(file.path(folder,'FAILED.txt'))) file.rename(file.path(folder,'FAILED.txt'),
          file.path(folder,paste0('PREVIOUS_',format(Sys.time(),'%Y%m%d_%H%M%S'),'_FAILED.txt')))
        results[[g]]<-ans
        run_status[[g]]<-data.table::data.table(protein=g,status=.pcfr_run_status(ans),
          detail=if(g=='CD101' && .pcfr_run_status(ans)=='assessment_complete_no_PTSD_credible_set')
            'Both PTSD fits converged with no 95% credible set; protein fit not run.' else '')
      }
      .pcfr_aggregate(results,out)
      .pcfr_write(data.table::rbindlist(run_status,fill=TRUE),file.path(out,'run_status.csv'))
      gc(verbose=FALSE)
    }
    cat('STEP 4/4: Final review summaries\n')
    status<-data.table::rbindlist(run_status,fill=TRUE)
    writeLines(c('Computation loop finished. Scientific interpretation requires review of coverage, convergence, LD consistency and all signal pairs.',
      paste('Completed regions:',sum(status$status %in% c('completed','assessment_complete_no_PTSD_credible_set')),'of',nrow(targets)),
      paste('Regions needing technical review:',sum(!status$status %in% c('completed','assessment_complete_no_PTSD_credible_set')))),file.path(out,'FINISHED.txt'))
    print(status)
    if(length(results)) print(data.table::rbindlist(lapply(results,function(x)x$overview),fill=TRUE))
    controls<-lapply(results,function(x)x$control);controls<-controls[vapply(controls,function(x)!is.null(x),logical(1))]
    if(length(controls) && any(!data.table::rbindlist(controls,fill=TRUE)$matches))
      cat('Chr15 control differences were detected. Review the control comparison; no released results were replaced.\n')
    cat('Output directory: ',out,'\n',sep='')
  },error=function(e) {
    fatal<<-conditionMessage(e);writeLines(fatal,file.path(out,'FAILED.txt'))
    cat('STOPPED: ',fatal,'\nOutput directory: ',out,'\n',sep='')
  })
  if(!is.null(fatal)) stop(fatal,call.=FALSE)
  invisible(out)
}

if(!isTRUE(getOption('ptsd.coloc.full_region.load_only',FALSE))) ptsd_coloc_full_region()
