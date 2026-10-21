# Frozen dense-region colocalization. No model fitting or automatic snapshot update.
(function() {
  P <- getOption('ptsd.current.paths');if(is.null(P))stop('Source 00_START.R first')
  input <- file.path(P$snapshot_dir,'coloc_full_region_20261020')
  read <- function(n) utils::read.csv(file.path(input,n),check.names=FALSE,stringsAsFactors=FALSE)
  su <- read('coloc_SuSiE_prior_grid.csv');ab <- read('coloc_ABF_prior_grid.csv')
  q <- read('summary_counts.csv');overview <- read('reference_prior_overview.csv')
  stopifnot(nrow(q)==10L,nrow(ab)==70L,nrow(su)==462L,all(is.finite(ab$PP.H4.abf)))
  mx <- stats::aggregate(su$PP.H4.abf,list(protein=su$protein,N_mode=su$N_mode,p12=su$p12),max)
  names(mx)[4] <- 'max_signal_pair_PP_H4'
  out <- file.path(P$report_dir,'Supplementary_Figure_S9');dir.create(out,recursive=TRUE,showWarnings=FALSE)
  utils::write.csv(mx,file.path(out,'SuSiE_maximum_pair_plot_data.csv'),row.names=FALSE)
  utils::write.csv(ab,file.path(out,'ABF_plot_data.csv'),row.names=FALSE)
  genes <- c('FES','FURIN','CGREF1','KHK','UBE2L6','CD40','AKT3','SIRPA','SNX18','CD101')
  cols <- c(NTOT='#2166AC',NEFF='#00897B',ABF='#C46B19')
  draw <- function() {
    graphics::par(mfrow=c(5,2),mar=c(2.7,3.8,2.3,1.1),oma=c(2.4,0,3.2,0),
      family='sans',bg='white',mgp=c(2.2,.65,0),cex=.93)
    for(g in genes) {
      n <- q$final_variants[match(g,q$protein)]
      graphics::plot(NA,xlim=c(1e-7,1e-4),ylim=c(0,1),log='x',xaxt='n',
        xlab='',ylab='Posterior H4',main=paste0(g,'  (',format(n,big.mark=','),' SNPs)'),
        cex.main=1.05,las=1)
      graphics::axis(1,at=c(1e-7,1e-6,1e-5,1e-4),labels=c('1e-7','1e-6','1e-5','1e-4'))
      graphics::abline(v=1e-5,lty=3,col='grey65');graphics::abline(h=.8,lty=2,col='grey75')
      d <- ab[ab$protein==g,];d<-d[order(d$p12),]
      graphics::lines(d$p12,d$PP.H4.abf,type='o',pch=17,col=cols['ABF'],lwd=1.3,cex=.6)
      for(mode in c('NTOT','NEFF')) {
        d<-mx[mx$protein==g & mx$N_mode==mode,];d<-d[order(d$p12),]
        if(nrow(d))graphics::lines(d$p12,d$max_signal_pair_PP_H4,type='o',
          pch=if(mode=='NTOT')19 else 1,lty=if(mode=='NTOT')1 else 2,
          col=cols[mode],lwd=1.4,cex=if(mode=='NTOT').65 else .95)
      }
      if(!g %in% unique(su$protein)) {
        msg<-if(g=='SIRPA') 'SuSiE: protein fit did not converge' else 'SuSiE: no PTSD credible set'
        graphics::mtext(msg,side=3,line=.2,cex=.65,col='#555555')
      }
    }
    graphics::mtext('Dense regional colocalization and prior sensitivity',side=3,outer=TRUE,line=1.8,font=2,cex=1.2)
    graphics::mtext('SuSiE NTOT: blue solid   |   SuSiE NEFF: green dashed   |   ABF: orange triangles',side=3,outer=TRUE,line=.35,cex=.83)
    graphics::mtext('Shared-variant prior p12 (p1 = p2 = 1e-4)',side=1,outer=TRUE,line=.5,cex=.95)
    graphics::mtext('SuSiE curves show the maximum across signal pairs. Missing curves are unavailable, not zero.',side=1,outer=TRUE,line=1.8,cex=.78)
  }
  formats <- getOption('ptsd.coloc.plot.formats',c('pdf','png','tiff'))
  stopifnot(length(formats)>0L,all(formats %in% c('pdf','png','tiff')))
  for(f in unique(formats)) {
    path<-file.path(out,paste0('SuppFig9_colocalization_current.',f))
    if(f=='pdf')grDevices::pdf(path,width=10,height=12,bg='white',useDingbats=FALSE)
    if(f=='png')grDevices::png(path,width=10,height=12,units='in',res=400,bg='white')
    if(f=='tiff')grDevices::tiff(path,width=10,height=12,units='in',res=400,compression='lzw',bg='white')
    tryCatch(draw(),finally=grDevices::dev.off())
  }
  caption<-paste('Supplementary Figure S9. Dense-region colocalization prior sensitivity for all ten proteins.',
    'Each panel uses the full matched, QC-passing sentinel +/-500-kb regional set and the same seven p12 values with p1=p2=1e-4.',
    'Blue and green curves show maximum SuSiE signal-pair H4 under NTOT and source-NEFF settings; orange triangles show single-causal-variant ABF H4.',
    'The vertical dotted line marks reference p12=1e-5; the horizontal dashed line marks H4=0.8.',
    'AKT3, CD40, SNX18 and CD101 have no PTSD 95% credible set under either N setting. SIRPA protein fitting did not converge, and its PTSD fits also had no credible set.',
    'CD101 protein fitting was explicitly not run after outcome assessment. Missing curves are not zeros or evidence against sharing.',
    'CD40 ABF support is method-dependent. CGREF1 and KHK remain prior-sensitive; the maximum need not involve the sentinel signal.',
    'CD40 and CGREF1 protein fits reached the L=10 credible-set cap. External LD mismatch and regional median N approximations limit resolution.',
    'All pairwise results are in Table S11c; coverage and fit diagnostics are in Table S11d.')
  writeLines(caption,file.path(out,'SuppFig9_caption.txt'),useBytes=TRUE)
  cat('Supplementary Figure S9 written from ten dense-region snapshots. No models rerun.\n')
})()
