# Direction-aware Figure 4. Plot only; no GWAS/MR/SMR analysis is rerun.
(function() {
  P <- getOption('ptsd.current.paths')
  output_parent <- getOption('ptsd.figure4.output_dir', if(is.null(P)) getwd() else P$report_dir)
  formats <- getOption('ptsd.figure4.formats',c('pdf','png','tiff'))
  stopifnot(length(formats)>0,all(formats %in% c('pdf','png','tiff')))
  if(is.null(P)) stop('Source 00_START.R first')
  d <- utils::read.csv(file.path(P$snapshot_dir,'sensitivity_review_20261019/eqtl_smr_direction.csv'),stringsAsFactors=FALSE,na.strings=c('NA',''))
  d$qualified <- as.logical(d$qualified)
  stopifnot(nrow(d)==20L, !anyDuplicated(paste(d$gene,d$tissue)))
  stopifnot(sum(d$qualified)==9L,
    sum(d$qualified & d$direction=='Opposite')==5L,
    sum(is.na(d$b_SMR))==6L)
  genes <- c('UBE2L6','SNX18','KHK','FES','CGREF1','SIRPA','FURIN','CD40','AKT3','CD101')
  d$x <- match(d$tissue,c('blood','brain'));d$y <- 11-match(d$gene,genes)
  stopifnot(all(is.finite(d$x)),all(is.finite(d$y)))
  d$signed_log10_P <- sign(d$b_SMR)*(-log10(d$p_SMR))
  palette <- grDevices::colorRampPalette(c('#2166AC','#67A9CF','#F7F7F7','#EF8A62','#B2182B'))(257)
  limit <- ceiling(max(abs(d$signed_log10_P),na.rm=TRUE))
  colors <- function(z) {
    ans <- rep('#E4E7EB',length(z));ok <- is.finite(z)
    ans[ok] <- palette[pmax(1,pmin(257,1+round((z[ok]+limit)/(2*limit)*256)))]
    ans
  }
  draw <- function() {
    graphics::par(omi=c(.4,.35,.4,.4),bg='white',family='sans')
    graphics::layout(matrix(c(1,2,3,3),2,byrow=TRUE),widths=c(7.1,1.0),heights=c(6.2,1.5))
    graphics::par(mar=c(.2,5.1,4.4,.6),xaxs='i',yaxs='i')
    graphics::plot.new();graphics::plot.window(xlim=c(.5,3.45),ylim=c(.5,10.5))
    graphics::rect(d$x-.48,d$y-.48,d$x+.48,d$y+.48,col=colors(d$signed_log10_P),border='white')
    graphics::axis(2,at=10:1,labels=genes,las=1,tick=FALSE,line=-.3,cex.axis=1.00,font=2)
    graphics::axis(3,at=1:3,labels=c('Blood\n(Westra)','Brain\n(BrainMeta v2)','Plasma protein\nMR sign'),tick=FALSE,line=.1,cex.axis=.92,font=2)
    graphics::mtext('Transcript association direction and plasma protein MR',side=3,line=3.0,adj=0,cex=1.08,font=2)
    tc <- ifelse(abs(d$signed_log10_P)>3.7,'white','#172B40')
    for(i in seq_len(nrow(d))) {
      if(!is.finite(d$b_SMR[i])) {graphics::text(d$x[i],d$y[i],'NA',col='#626B75',cex=.9);next}
      marker <- if(d$qualified[i]) ' *' else if(d$statistical_result[i]=='HEIDI rejected') ' # ' else ''
      graphics::text(d$x[i],d$y[i]+.09,paste0(sprintf('%+.4f',d$b_SMR[i]),marker),col=tc[i],cex=.94,font=2)
      label <- if(d$qualified[i] && d$direction[i]=='Opposite') 'opposite' else if(d$statistical_result[i]=='HEIDI rejected') 'HEIDI rejected' else ''
      if(nzchar(label)) graphics::text(d$x[i],d$y[i]-.21,label,col=tc[i],cex=.78)
    }
    for(g in genes) {
      b <- d$plasma_MR_beta[match(g,d$gene)];y <- 11-match(g,genes)
      graphics::rect(2.68,y-.48,3.43,y+.48,col='#F8F9FB',border='white')
      graphics::text(3.05,y,if(b>0) '+' else '-',col=if(b>0) '#B2182B' else '#2166AC',cex=1.6,font=2)
    }
    graphics::par(mar=c(.2,.2,4.4,3.1),xaxs='i',yaxs='i')
    graphics::plot.new();graphics::plot.window(xlim=c(0,1),ylim=c(-limit,limit))
    ed <- seq(-limit,limit,length.out=258)
    graphics::rect(.05,ed[-258],.6,ed[-1],col=palette,border=NA)
    graphics::axis(4,at=seq(-limit,limit,by=1),las=1,cex.axis=.76,tck=-.08,col='#CBD5E1')
    graphics::mtext('Signed\n-log10(P)',side=3,line=.7,cex=.83)
    graphics::par(mar=c(.4,.7,.2,.6),xaxs='i',yaxs='i')
    graphics::plot.new();graphics::plot.window(xlim=c(0,1),ylim=c(0,1))
    lines <- c('*  SMR FDR < 0.05 and HEIDI P > 0.05 (>=3 HEIDI SNPs).',
      '#  SMR FDR < 0.05 and HEIDI P <= 0.05 (possible linkage).',
      'Red: positive transcript effect. Blue: negative transcript effect. Numbers: SMR beta.',
      '"opposite": qualified transcript effect has the opposite sign to plasma protein MR.',
      'NA: SMR unavailable. Compare direction across layers; effect scales differ.')
    graphics::text(.01,seq(.90,.10,length.out=5),lines,adj=0,cex=.79,col='#273444')
  }
  dir.create(output_parent,recursive=TRUE,showWarnings=FALSE)
  out <- tempfile(paste0('Figure4_',format(Sys.time(),'%Y%m%d_%H%M%S'),'_'),tmpdir=output_parent)
  stopifnot(dir.create(out))
  for(f in unique(formats)) {
    fn <- file.path(out,paste0('Figure4_eQTL_SMR_current.',f))
    if(f=='pdf') grDevices::pdf(fn,width=9.0,height=8.5,useDingbats=FALSE,bg='white')
    if(f=='png') grDevices::png(fn,width=9.0,height=8.5,units='in',res=600,bg='white',pointsize=12)
    if(f=='tiff') grDevices::tiff(fn,width=9.0,height=8.5,units='in',res=600,compression='lzw',bg='white',pointsize=12)
    tryCatch(draw(),finally=grDevices::dev.off())
  }
  utils::write.csv(d[setdiff(names(d),c('x','y'))],file.path(out,'Figure4_plot_data.csv'),row.names=FALSE,na='NA')
  caption <- paste('Figure 4. Direction-aware brain and blood eQTL-SMR results for PTSD-prioritized genes.',
    'Fill represents sign(beta_SMR) multiplied by -log10(P_SMR); red is positive and blue is negative.',
    'Numbers show transcript-level SMR coefficients. The right column shows only the sign of plasma pQTL-MR from Table S3a, on a different exposure scale.',
    'An asterisk marks within-tissue SMR FDR < 0.05 plus HEIDI P > 0.05 with at least three HEIDI SNPs; a hash marks FDR-significant SMR with HEIDI rejection.',
    'Opposite labels identify qualified transcript associations with a sign opposite to plasma protein MR: brain and blood KHK, brain SNX18 and UBE2L6, and blood FURIN.',
    'Brain FURIN is concordant, yielding mixed tissue directions. Brain FES is HEIDI-rejected; blood FES qualifies and is concordant.',
    'Gray NA cells indicate unavailable SMR. Statistical qualification and direction concordance are distinct; HEIDI non-rejection alone does not establish colocalization.',
    'Source: Additional file 3, Tables S14 and S3a. No association estimates were recalculated.')
  writeLines(caption,file.path(out,'Figure4_caption.txt'),useBytes=TRUE)
  options(ptsd.figure4.last_output=normalizePath(out,winslash='/',mustWork=TRUE))
  cat('Figure 4 written: 9 qualified tissue-gene associations; 5 opposite directions across 4 genes. Plot only.\nOutput: ',out,'\n',sep='')
  invisible(out)
})()
