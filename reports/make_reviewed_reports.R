# Generate revised result tables and draft figures from immutable reviewed snapshots.
(function() {
  P<-getOption("ptsd.current.paths");if(is.null(P)) stop("Source 00_START.R first")
  dir.create(P$report_dir,recursive=TRUE,showWarnings=FALSE)
  read <- function(f) read.csv(file.path(P$snapshot_dir,f),check.names=FALSE,stringsAsFactors=FALSE)
  save <- function(d,f) write.csv(d,file.path(P$report_dir,f),row.names=FALSE,na="NA")
  for(n in c("eqtl_smr_direction.csv","reverse_mr_review.csv","evidence_matrix.csv"))
    save(read(file.path("sensitivity_review_20261019",n)),paste0("reviewed_",n))
  factor<-read("sem_expansion_20261019/primary/factor_MR_expanded.csv")
  factor$effect_allele<-factor$exposure_EA;factor$other_allele<-factor$exposure_OA
  save(factor,"S26_factor_MR_corrected_allele_labels.csv")
  # Historical chr15 files at the snapshot root are retained only as controls.
  for(n in c("coloc_ABF_prior_grid.csv","coloc_SuSiE_prior_grid.csv","fit_status.csv",
             "summary_counts.csv","reference_prior_overview.csv","reviewed_summary.csv"))
    save(read(file.path("coloc_full_region_20261020",n)),paste0("dense_region_",n))
  m<-read("mediation_results.csv");save(m,"CRP_indirect_effects_both_windows.csv")
  n<-read("network_results.csv");save(n,"exploratory_network_six_pairs.csv")
  save(read("CRP_to_PTSD_results.csv"),"CRP_to_PTSD_both_windows.csv")
  # Two-step indirect estimates; FURIN unavailable, not zero.
  pdf(file.path(P$report_dir,"SuppFig6_CRP_indirect_effects_current.pdf"),width=10,height=5.5)
  par(mfrow=c(1,2),mar=c(5,6,4,1),oma=c(2,0,0,0))
  eligible <- m[m$status=="testable" & m$analysis %in% c("cis500kb","cis100kb"),]
  common_lim <- range(c(eligible$indirect_lower,eligible$indirect_upper,0),finite=TRUE)
  for(a in c("cis500kb","cis100kb")) {
    d<-m[m$analysis==a&m$status=="testable",];d<-d[order(d$protein),]
    lim<-common_lim
    plot(d$indirect_beta,seq_len(nrow(d)),xlim=lim,ylim=c(.5,nrow(d)+.5),yaxt="n",pch=19,
         xlab="Exploratory indirect estimate (PTSD log-odds scale)",ylab="",main=if(a=="cis500kb") "A  CRP cis +/-500 kb" else "B  CRP cis +/-100 kb")
    abline(v=0,col="grey65",lty=2);axis(2,at=seq_len(nrow(d)),labels=d$protein,las=1)
    segments(d$indirect_lower,seq_len(nrow(d)),d$indirect_upper,seq_len(nrow(d)))
  }
  mtext("95% delta-method intervals; zero covariance assumed. FURIN not evaluable (CRP EAF unavailable).",side=1,outer=TRUE,cex=.8)
  dev.off()
  pdf(file.path(P$report_dir,"Exploratory_network_current.pdf"),width=8,height=5.8)
  par(mar=c(7,10,4,2));y<-seq_len(nrow(n));lim<-range(c(n$lower,n$upper,0))
  plot(n$beta,y,xlim=lim,ylim=c(.5,nrow(n)+.5),yaxt="n",pch=19,xlab="Exploratory protein-to-protein Wald estimate",ylab="",main="Six selected protein pairs")
  abline(v=0,col="grey65",lty=2);axis(2,at=y,labels=paste(n$upstream,"->",n$downstream),las=1)
  segments(n$lower,y,n$upper,y);mtext("No pair passes BH-FDR < 0.05; sample overlap unresolved.",side=3,cex=.8)
  mtext("Intervals are unadjusted 95% confidence intervals.",side=1,line=5.2,cex=.8)
  dev.off()
  source(file.path(P$package_root,"reports/make_coloc_reports.R"),local=new.env(parent=globalenv()))
  cat("Reviewed component tables and plots exported; full Word/XLSX documents are not rewritten.\n")
})()
