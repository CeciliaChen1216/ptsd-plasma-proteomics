# Complete the CD101 colocalization assessment from the reviewed run's cache.
# Revision 2026.10.20.2. Source this file; no original analysis is overwritten.
# Default previous run: coloc_full_region_20261019_163303_c9827f33c6b.
# Outcome fits run first. If neither yields a 95% credible set, protein fitting
# is unnecessary for coloc.susie and is explicitly recorded as not run.
# All fitted models retain the original SNPs, LD, N conventions and parameters.
# R.utils applies a 15-minute time limit per fit, subject to R interrupt support.
# This is a cache completion helper, not a replacement for the main pipeline.

.pcfc_no_outcome_cs <- function(fits) {
  length(fits)==2L && all(vapply(fits,function(a)
    !is.null(a$fit) && length(a$fit$sets$cs)==0L,logical(1)))
}

.pcfc_pair_status <- function(protein,outcome) {
  if(is.null(outcome$fit)) 'fit_error_or_nonconvergence' else
    if(!length(outcome$fit$sets$cs)) 'no_PTSD_credible_set' else
    if(is.null(protein$fit)) 'fit_error_or_nonconvergence' else
    if(!length(protein$fit$sets$cs)) 'no_protein_credible_set' else ''
}

.pcfc_read_abf_cache <- function(path,prior_grid) {
  # With automatic type inference, an entirely blank error column can become
  # logical NA; nzchar(NA) is TRUE. Preserve the original character column.
  ab <- data.table::fread(path,colClasses=list(character='error'))
  required <- c('protein','p12','error','PP.H4.abf')
  if(!all(required %in% names(ab))) stop('The CD101 ABF cache lacks required columns')
  if(nrow(ab)!=length(prior_grid) || !setequal(ab$p12,prior_grid) ||
     !all(ab$protein=='CD101') || !all(is.finite(ab$PP.H4.abf)) ||
     anyNA(ab$error) || any(nzchar(trimws(ab$error))))
    stop('The CD101 ABF cache is incomplete or contains errors')
  ab
}

.pcfc_run <- function() {
  previous <- getOption('ptsd.coloc.completion.previous_dir',
    '')
  pipeline <- getOption('ptsd.coloc.completion.pipeline',
    file.path(getOption('ptsd.current.root',''),
      'provenance/executed_scripts/PTSD_coloc_full_region_rerun.R.txt'))
  if(!nzchar(previous)) stop('Set options(ptsd.coloc.completion.previous_dir=...) to the original run directory first')
  expected_md5 <- 'bb93c219b308f1236776fcc84a23d537'
  paths <- c(state=file.path(previous,'RUN_STATE.rds'),pipeline=pipeline,
    target=file.path(previous,'target_manifest.csv'))
  if(any(!file.exists(paths))) stop('Missing original run/cache files: ',
    paste(paths[!file.exists(paths)],collapse='; '))
  state <- readRDS(paths[['state']])
  if(!identical(unname(tools::md5sum(pipeline)),expected_md5) ||
     !identical(state$script_md5,expected_md5))
    stop('Original pipeline differs from the reviewed 2026.10.19.1 script. Preserve the original file and request review.')
  cfg <- state$config
  need <- c('data.table','R.utils','coloc','susieR')
  missing <- need[!vapply(need,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing)) stop('Missing packages: ',paste(missing,collapse=', '))
  versions <- data.frame(package=c('R',need),version=c(as.character(getRversion()),
    vapply(need,function(p)as.character(utils::packageVersion(p)),character(1))))
  saved_versions <- as.data.frame(state$package_versions)
  if(!identical(versions$package,saved_versions$package) ||
     !identical(versions$version,saved_versions$version))
    stop('R/package versions changed. Use the original versions or review this as a separate sensitivity analysis.')
  suppressPackageStartupMessages(library(data.table))
  old_load <- getOption('ptsd.coloc.full_region.load_only')
  options(ptsd.coloc.full_region.load_only=TRUE)
  on.exit(options(ptsd.coloc.full_region.load_only=old_load),add=TRUE)
  engine <- new.env(parent=globalenv())
  sys.source(pipeline,envir=engine)
  targets <- data.table::fread(paths[['target']])
  if(!setequal(targets$protein,c('FES','FURIN','AKT3','CD40','CGREF1','KHK',
    'SIRPA','SNX18','UBE2L6','CD101'))) stop('Unexpected target manifest')
  other <- setdiff(targets$protein,'CD101')
  checkpoints <- file.path(previous,'regions',other,'REGION_COMPLETE.rds')
  if(any(!file.exists(checkpoints)))
    stop('An earlier protein checkpoint is missing: ',paste(other[!file.exists(checkpoints)],collapse=', '))
  out <- tempfile(paste0('coloc_CD101_completion_',format(Sys.time(),'%Y%m%d_%H%M%S'),'_'),
    tmpdir=cfg$output_parent)
  if(!dir.create(out,recursive=TRUE)) stop('Cannot create completion output')
  out <- normalizePath(out,winslash='/',mustWork=TRUE)
  cat('Completion output directory: ',out,'\n',sep='')
  log <- file(file.path(out,'RUN.log'),'wt')
  sink(log,split=TRUE);sink(log,type='message')
  on.exit({
    capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
    sink(type='message');sink();close(log)
    tryCatch(engine$.pcfr_sendback(out),error=function(e)
      message('Upload the SEND_BACK folder if archive creation failed: ',conditionMessage(e)))
    if(file.exists(file.path(out,'FAILED.txt')))
      cat('STATUS: STOPPED. The archive is a diagnostic bundle; inspect FAILED.txt.\n')
    else if(!file.exists(file.path(out,'FINISHED.txt')))
      cat('STATUS: INCOMPLETE. No completion marker was written.\n')
    else cat('STATUS: ASSESSMENT COMPLETE. Review the saved fit and colocalization statuses.\n')
  },add=TRUE)
  tryCatch({
    cat('STEP 1/3: Reuse completed regions and verify CD101 caches\n')
    send <- file.path(out,'SEND_BACK');dir.create(send)
    if(dir.exists(file.path(previous,'SEND_BACK'))) {
      old_small <- list.files(file.path(previous,'SEND_BACK'),full.names=TRUE)
      old_small <- old_small[!file.info(old_small)$isdir &
        grepl('^(FES|FURIN|AKT3|CD40|CGREF1|KHK|SIRPA|SNX18|UBE2L6)__',basename(old_small))]
      if(length(old_small)) file.copy(old_small,file.path(send,basename(old_small)))
    }
    for(f in c('input_manifest.csv','package_versions.csv','target_manifest.csv',
      'PTSD_source_header.txt','config_used.txt','reviewed_chr15_control_comparison.csv'))
      if(file.exists(file.path(previous,f))) file.copy(file.path(previous,f),file.path(out,f))
    file.copy(file.path(previous,'RUN.log'),file.path(send,'PREVIOUS_RUN.log'))
    file.copy(pipeline,file.path(send,'BASE_PIPELINE.R.txt'))
    source_files <- vapply(sys.frames(),function(x) if(is.character(x$ofile)&&length(x$ofile)==1L) x$ofile else '',character(1))
    source_files <- source_files[nzchar(source_files)]
    if(length(source_files)) file.copy(tail(source_files,1),file.path(out,'script_snapshot.R.txt'))
    results <- stats::setNames(lapply(checkpoints,readRDS),other)
    target <- targets[protein=='CD101']
    old_region <- file.path(previous,'regions','CD101')
    folder <- file.path(out,'regions','CD101');dir.create(folder,recursive=TRUE)
    old_checkpoint <- file.path(old_region,'REGION_COMPLETE.rds')
    original_status <- data.table::fread(file.path(previous,'run_status.csv'))
    if(file.exists(old_checkpoint)) {
      cat('CD101 has subsequently completed in the original directory; reusing its completed checkpoint.\n')
      results[['CD101']] <- readRDS(old_checkpoint)
      small <- list.files(old_region,pattern='[.]csv$|credible_sets[.]txt$|LD_consistency[.]txt$',full.names=TRUE)
      small <- small[!grepl('analysis_input|preLD_audit|frequency_audit|reference_QC',basename(small))]
      file.copy(small,file.path(folder,basename(small)))
      cd_status <- original_status[protein=='CD101']
      if(!nrow(cd_status)) cd_status <- data.table::data.table(protein='CD101',
        status='completed_checkpoint_reused',detail='Review the saved fit_status and colocalization status.')
    } else {
      caches <- c(analysis=file.path(old_region,'analysis_input.csv'),LD=file.path(old_region,'LD.rds'),
        counts=file.path(old_region,'summary_counts.csv'),ABF=file.path(old_region,'coloc_ABF_prior_grid.csv'))
      if(any(!file.exists(caches))) stop('A CD101 cache is missing: ',paste(caches[!file.exists(caches)],collapse='; '))
      cache_manifest <- data.table::data.table(item=names(caches),path=unname(caches),
        md5=unname(tools::md5sum(caches)))
      engine$.pcfr_write(cache_manifest,file.path(out,'cache_manifest.csv'))
      d <- data.table::fread(caches[['analysis']]);R <- readRDS(caches[['LD']])
      if(nrow(d)!=1837L || anyDuplicated(d$SNP) || !all(d$reason=='keep') ||
         !identical(rownames(R),d$SNP) || !identical(colnames(R),d$SNP) ||
         !all(is.finite(R)) || max(abs(R-t(R)))>1e-8 || max(abs(diag(R)-1))>1e-8)
        stop('CD101 cache count, QC flags or LD order differ from the reviewed 1837-variant dataset')
      if(!all(d$CHR==1L & abs(d$BP-117544331L)<=500000L)) stop('CD101 cache coordinates differ')
      ds <- engine$.pcfr_datasets(d)
      if(ds$x$N!=33796) stop('Unexpected CD101 protein sample size')
      catalog <- data.table::as.data.table(readRDS(cfg$pqtl_rds))
      engine$.pcfr_sentinel(d,target,catalog,folder)
      counts <- data.table::fread(caches[['counts']])
      ab <- .pcfc_read_abf_cache(caches[['ABF']],cfg$p12_grid)
      engine$.pcfr_write(ab,file.path(folder,'coloc_ABF_prior_grid.csv'))
      for(f in c('CD101_LD_consistency.txt','preLD_reason_counts.csv','frequency_reason_counts.csv'))
        if(file.exists(file.path(old_region,f))) file.copy(file.path(old_region,f),file.path(folder,f))
      engine$.pcfr_write(counts,file.path(folder,'summary_counts.csv'))
      original_fit <- engine$.pcfr_fit
      bounded_fit <- function(z,n,label) {
        begin <- proc.time()[['elapsed']]
        ans <- tryCatch(R.utils::withTimeout(original_fit(z,R,n,label,cfg,folder),
          timeout=900,onTimeout='error'),error=function(e) {
            engine$.pcfr_save(e,file.path(folder,paste0(label,'_fit.rds')))
            list(fit=NULL,status=data.table::data.table(fit=label,n=n,converged=FALSE,
              nCS=NA_integer_,error=conditionMessage(e),warnings=''))
          })
        ans$status[,`:=`(execution='executed',elapsed_seconds=proc.time()[['elapsed']]-begin)]
        engine$.pcfr_write(ans$status,file.path(folder,paste0(label,'_status.csv')))
        ans
      }
      cat('STEP 2/3: Fit CD101 PTSD outcomes first, retaining all 1837 SNPs\n')
      fy <- lapply(c('NTOT','NEFF'),function(mode)
        bounded_fit(ds$zy,median(d[[mode]]),paste0('CD101_PTSD_',mode)))
      names(fy) <- c('NTOT','NEFF')
      # A missing PTSD credible set prevents coloc.susie regardless of the
      # protein fit. This shortcut changes computational order, not SNP/model QC.
      no_outcome_cs <- .pcfc_no_outcome_cs(fy)
      if(no_outcome_cs) {
        cat('Both PTSD fits converged without a 95% credible set. Protein refitting is not required to establish that SuSiE colocalization is unavailable.\n')
        fx <- list(fit=NULL,status=data.table::data.table(fit='CD101_protein',n=ds$x$N,
          converged=NA,nCS=NA_integer_,error='',warnings='',
          execution='not_run_no_PTSD_credible_set',elapsed_seconds=0))
        engine$.pcfr_write(fx$status,file.path(folder,'CD101_protein_status.csv'))
      } else fx <- bounded_fit(ds$zx,ds$x$N,'CD101_protein')
      fits <- data.table::rbindlist(c(list(fx$status),lapply(fy,function(a)a$status)),fill=TRUE)
      pairs <- list();statuses <- list()
      for(mode in c('NTOT','NEFF')) for(pr in cfg$p12_grid) {
        key <- paste(mode,pr,sep='_');err<-'';npairs<-0L
        why <- .pcfc_pair_status(fx,fy[[mode]])
        if(nzchar(why)) status<-why else {
          cs<-tryCatch(coloc::coloc.susie(fx$fit,fy[[mode]]$fit,
            p1=cfg$p1,p2=cfg$p2,p12=pr),error=function(e)e)
          engine$.pcfr_save(cs,file.path(folder,paste0('coloc_',key,'.rds')))
          if(inherits(cs,'error')) {status<-'coloc_error';err<-conditionMessage(cs)} else
            if(is.null(cs$summary)||!nrow(cs$summary)) status<-'no_signal_pairs_returned' else {
              status<-'estimated';npairs<-nrow(cs$summary)
              pairs[[key]]<-cbind(data.table::data.table(protein='CD101',N_mode=mode,p12=pr),
                data.table::as.data.table(cs$summary))
            }
        }
        statuses[[key]]<-data.table::data.table(protein='CD101',N_mode=mode,p12=pr,
          status=status,n_pairs=npairs,error=err)
      }
      su<-if(length(pairs)) data.table::rbindlist(pairs,fill=TRUE) else NULL
      ps<-data.table::rbindlist(statuses,fill=TRUE)
      overview<-data.table::rbindlist(lapply(c('NTOT','NEFF'),function(mode) {
        z<-if(is.null(su)) NULL else su[N_mode==mode & abs(p12-cfg$p12_reference)<1e-15]
        v<-if(is.null(z)||!nrow(z)) NULL else z[which.max(PP.H4.abf)]
        ss<-ps[N_mode==mode & abs(p12-cfg$p12_reference)<1e-15]
        data.table::data.table(protein='CD101',N_mode=mode,n_snps=nrow(d),p12=cfg$p12_reference,
          nCS_protein=fx$status$nCS,nCS_PTSD=fy[[mode]]$status$nCS,
          status=ss$status[1],n_pairs=ss$n_pairs[1],
          maximum_pair_PP_H4=if(is.null(v)) NA_real_ else v$PP.H4.abf,
          same_pair_PP_H3=if(is.null(v)) NA_real_ else v$PP.H3.abf,
          hit1=if(is.null(v)) NA_character_ else as.character(v$hit1),
          hit2=if(is.null(v)) NA_character_ else as.character(v$hit2))
      }),fill=TRUE)
      ans<-list(counts=counts,fit_status=fits,abf=ab,susie=su,
        pair_status=ps,overview=overview,control=NULL)
      results[['CD101']]<-ans
      for(item in c('fit_status','pair_status','susie')) {
        filename<-c(fit_status='fit_status.csv',pair_status='coloc_status.csv',susie='coloc_SuSiE_prior_grid.csv')[[item]]
        engine$.pcfr_write(ans[[item]],file.path(folder,filename))
      }
      engine$.pcfr_save(ans,file.path(folder,'REGION_COMPLETE.rds'))
      has_error<-any(ps$status %in% c('fit_error_or_nonconvergence','coloc_error','no_signal_pairs_returned'))
      cd_status<-data.table::data.table(protein='CD101',status=if(has_error) 'completed_with_fit_or_coloc_error' else
        if(no_outcome_cs) 'assessment_complete_no_PTSD_credible_set' else 'completed',
        detail=if(no_outcome_cs) 'Both PTSD fits converged; protein refit skipped because neither PTSD fit yielded a 95% credible set.' else '')
    }
    cat('STEP 3/3: Combine all ten assessments, retaining the prior SIRPA nonconvergence\n')
    engine$.pcfr_aggregate(results,out)
    status<-data.table::rbindlist(list(original_status[protein!='CD101'],cd_status),fill=TRUE)
    engine$.pcfr_write(status,file.path(out,'run_status.csv'))
    note<-c('Cache completion 2026.10.20.2; original inputs and results are read-only.',
      paste('Previous run:',previous),
      'The original CD101 protein fit has no saved completion in the reviewed bundle; no cause of interruption is assumed.',
      'CD101 PTSD NTOT/NEFF fits run before the protein fit. If both converge without a 95% credible set, protein refitting is explicitly skipped.',
      'This computational shortcut does not alter the variant set, LD, parameters or colocalization decision; no PTSD credible set means coloc.susie is unavailable.',
      'Each newly executed fit has a 900-second R time limit, subject to R interrupt support. A time limit or failed fit is not negative biological evidence.',
      'No diagnostic-driven SNP removal, LD ridge adjustment, lowered credible-set coverage or prior selection is performed.',
      'SIRPA retains its original nonconvergence status. Single-causal-variant ABF and multi-signal SuSiE results must be reported separately.')
    writeLines(note,file.path(out,'parameters.txt'));writeLines(note,file.path(send,'COMPLETION_NOTE.txt'))
    if(file.exists(file.path(out,'cache_manifest.csv')))
      file.copy(file.path(out,'cache_manifest.csv'),file.path(send,'cache_manifest.csv'))
    writeLines(c('All ten colocalization assessments are represented; this does not imply all fine-mapping models converged.',
      'Review fit_status.csv, coloc_status.csv and the computational shortcut noted in parameters.txt.'),
      file.path(out,'FINISHED.txt'))
    print(status);print(data.table::rbindlist(lapply(results,function(a)a$overview),fill=TRUE))
    cat('Completion output directory: ',out,'\n',sep='')
    invisible(out)
  },error=function(e) {
    writeLines(conditionMessage(e),file.path(out,'FAILED.txt'))
    cat('STOPPED: ',conditionMessage(e),'\nOutput directory: ',out,'\n',sep='')
    stop(conditionMessage(e),call.=FALSE)
  })
}

if(!isTRUE(getOption('ptsd.coloc.cd101.finish.load_only',FALSE))) .pcfc_run()
