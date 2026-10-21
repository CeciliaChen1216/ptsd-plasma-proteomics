# Current Supplementary Figure S2. Reviewed categories, not an evidence score.
(function() {
  P <- getOption('ptsd.current.paths');if(is.null(P)) stop('Source 00_START.R first')
  d <- utils::read.csv(file.path(P$snapshot_dir,'sensitivity_review_20261019/evidence_matrix.csv'),stringsAsFactors=FALSE)
  stopifnot(nrow(d)==80L,!anyDuplicated(paste(d$protein,d$domain)))
  genes <- c('AKT3','CGREF1','FES','KHK','SNX18','CD40','FURIN','UBE2L6','SIRPA','CD101')
  domains <- unique(d$domain)
  labels <- c('Primary\nMR','Steiger\ndirectionality','MR-PRESSO\nassessment','Reverse-MR\nassessment','Multi-SNP\nIVW','eQTL-SMR\n+ direction','Colocalization','Brain pQTL\n/ PWAS')
  d$x <- match(d$domain,domains);d$y <- 11-match(d$protein,genes)
  pal <- c(S='#2C7FB8',C='#D95F0E',I='#BDBDBD',N='#F0F0F0')
  out <- file.path(P$report_dir,'Supplementary_Figure_S2');dir.create(out,recursive=TRUE,showWarnings=FALSE)
  formats <- getOption('ptsd.suppfig2.formats',c('pdf','png','tiff'))
  stopifnot(length(formats)>0,all(formats %in% c('pdf','png','tiff')))
  draw <- function() {
    graphics::par(omi=c(.20,.15,.25,.20),bg='white',family='sans')
    graphics::layout(matrix(1:2,2),heights=c(5.3,.6))
    graphics::par(mar=c(.1,5.0,4.5,.5),xaxs='i',yaxs='i')
    graphics::plot.new();graphics::plot.window(xlim=c(.5,8.5),ylim=c(.5,10.5))
    graphics::rect(d$x-.485,d$y-.48,d$x+.485,d$y+.48,col=unname(pal[d$state]),border='white')
    graphics::axis(2,at=10:1,labels=genes,tick=FALSE,las=1,line=-.3,font=2,cex.axis=1.05)
    graphics::axis(3,at=1:8,labels=labels,tick=FALSE,line=.05,cex.axis=.89)
    graphics::mtext('Evidence matrix for ten PTSD-prioritized proteins',side=3,line=3.2,adj=0,cex=1.15,font=2)
    a <- which(!is.na(d$annotation)&nzchar(d$annotation))
    graphics::text(d$x[a],d$y[a],d$annotation[a],col=ifelse(d$state[a] %in% c('S','C'),'white','#283442'),cex=.83)
    graphics::par(mar=c(.1,.5,.1,.5));graphics::plot.new();graphics::plot.window(xlim=c(0,1),ylim=c(0,1))
    graphics::legend('center',legend=c('Supportive','Cautionary','Inconclusive','Not evaluable'),fill=pal,border=NA,horiz=TRUE,bty='n',cex=1.02,x.intersp=.7)
  }
  for(f in unique(formats)) {
    path <- file.path(out,paste0('SuppFig2_evidence_matrix_current.',f))
    if(f=='pdf')grDevices::pdf(path,width=10.4,height=6.35,bg='white',useDingbats=FALSE)
    if(f=='png')grDevices::png(path,width=10.4,height=6.35,units='in',res=600,bg='white',pointsize=12)
    if(f=='tiff')grDevices::tiff(path,width=10.4,height=6.35,units='in',res=600,compression='lzw',bg='white',pointsize=12)
    tryCatch(draw(),finally=grDevices::dev.off())
  }
  utils::write.csv(d[setdiff(names(d),c('x','y'))],file.path(out,'SuppFig2_plot_data.csv'),row.names=FALSE)
  caption <- paste('Supplementary Figure S2. Evidence matrix for ten PTSD-prioritized proteins.',
    'The eight domains are interpreted separately and are not summed.',
    'Reverse-MR cells are cautionary for all six proteins with IVW P < 0.05: CD40 has concordant IVW and weighted-median (WM) results; CGREF1, FES, KHK, SIRPA and UBE2L6 have IVW-only significance.',
    'Neither estimator is significant for the other four proteins; these results are inconclusive.',
    'Multi-instrument IVW was tested for SIRPA (P = 0.413, inconclusive) and CD101 (P = 0.221, cautionary because of the opposite estimate and heterogeneity).',
    'The eQTL-SMR domain combines statistical qualification (within-tissue FDR < 0.05, HEIDI P > 0.05, >=3 HEIDI SNPs) with sign comparison to plasma protein MR.',
    'KHK, SNX18 and UBE2L6 are cautionary for opposite qualified effects; FURIN is cautionary for mixed tissue directions. CD40, SIRPA and blood FES have concordant qualified results; brain FES is HEIDI-rejected.',
    'AKT3 and CD101 lack eligible transcript results; CGREF1 is tested but inconclusive.',
    'Colocalization uses uniform dense +/-500-kb regions. FES is supportive at p12=1e-5 (0.996*) but prior-sensitive. CGREF1 (0.773) and KHK (0.471) are inconclusive; FURIN and UBE2L6 favor distinct signals. CD40 is cautionary because ABF support (0.930) is not assessable by SuSiE.',
    'No CS denotes no PTSD 95% credible set for AKT3, SNX18 and CD101. SIRPA is not evaluable because the protein fit did not converge; both PTSD fits also had no credible set. CD101 protein fitting was not run. These states are not evidence of absence.',
    'Multi-instrument analyses used candidate-window credible-set lead variants. Brain-pQTL/PWAS states use their domain-specific definitions.',
    'Sources: Tables S1, S3a, S7b, S8, S11, S14 and S20. No absence-of-bias claim is inferred from a non-significant diagnostic.')
  writeLines(caption,file.path(out,'SuppFig2_caption.txt'),useBytes=TRUE)
  options(ptsd.suppfig2.last_output=normalizePath(out,winslash='/',mustWork=TRUE))
  cat('Supplementary Figure S2 written: IVW-only reverse signals and transcript direction conflicts are explicit. Plot only.\n')
  invisible(out)
})()
