# Figure 2: expanded, reviewed primary Genomic SEM results (2026-10-19).
# Plot only. No GWAS, LDSC, MR, or Genomic SEM model is rerun.
# No contributed R packages are required.
# Run after 00_START.R: ptsd_reports(figure2=TRUE)
# Optional overrides:
# options(ptsd.figure2.input="path/to/primary/factor_MR_expanded.csv")
# options(ptsd.figure2.output_dir="/path/to/figure-output")
# options(ptsd.figure2.formats=c("pdf","png","tiff"))

(function() {
  P <- getOption("ptsd.current.paths")
  if(is.null(P)) stop("Source 00_START.R first")
  root <- P$project_root
  input <- getOption("ptsd.figure2.input",file.path(P$snapshot_dir,"factor_MR.csv"))
  if(!file.exists(input)) stop("Missing reviewed primary snapshot: ",input)
  d <- utils::read.csv(input, stringsAsFactors=FALSE, check.names=FALSE)
  required <- c("protein","outcome_snp","mr_beta","mr_se","mr_p","mr_fdr","Q_pval","Q_fdr")
  if (!all(required %in% names(d))) stop("Input lacks required expanded factor-MR columns.")
  if (anyDuplicated(d$protein)) stop("Duplicate protein rows.")
  for (nm in c("mr_beta","mr_se","mr_p","mr_fdr","Q_pval","Q_fdr")) {
    if (!is.numeric(d[[nm]]) || any(!is.finite(d[[nm]]))) stop("Invalid values in ",nm)
  }
  if (any(d$mr_p<=0|d$mr_p>1|d$Q_pval<=0|d$Q_pval>1|d$mr_se<=0)) stop("Invalid P or SE values.")
  if (max(abs(d$mr_fdr-stats::p.adjust(d$mr_p,"BH")))>1e-10 ||
      max(abs(d$Q_fdr-stats::p.adjust(d$Q_pval,"BH")))>1e-10) stop("Stored BH FDRs do not match the full protein family.")
  candidates <- c("AKT3","CD40","CGREF1","FES","FURIN","KHK","SNX18","UBE2L6","CD101")
  counts <- c(proteins=nrow(d),unique_SNPs=length(unique(d$outcome_snp)),
    factor_MR_FDR_005=sum(d$mr_fdr<.05),Q_nominal_005=sum(d$Q_pval<.05),Q_FDR_005=sum(d$Q_fdr<.05),
    candidates=sum(d$protein %in% candidates),candidate_MR_FDR_005=sum(d$protein %in% candidates & d$mr_fdr<.05))
  expected <- c(1162L,1152L,51L,87L,0L,9L,8L)
  if (!identical(as.integer(counts),expected)) {
    print(counts)
    stop("Input does not match the reviewed primary analysis. Select primary/factor_MR_expanded.csv, not NEFF or an old snapshot.")
  }
  if (d$mr_fdr[d$protein=="CD101"]<.05 || "SIRPA" %in% d$protein) stop("Unexpected candidate coverage/status.")
  d$MR_minus_log10_P <- -log10(d$mr_p)
  d$Q_minus_log10_P <- -log10(d$Q_pval)
  d$MR_FDR_significant <- d$mr_fdr<.05
  d$PTSD_prioritized <- d$protein %in% candidates
  cd <- d[match(candidates,d$protein),,drop=FALSE]
  nominal <- -log10(.05)
  bonferroni <- -log10(.05/nrow(d))

  # The expanded SNP step reuses this recorded primary trait model unchanged.
  # Source: reviewed constrained_model_full.txt, 2026-10-07 primary LDSC run.
  model <- data.frame(trait=c("PTSD","MDD","Anxiety"),
    standardized_loading=c(.9260993,1.0039810,.8973164),
    standardized_residual=c(.1423401,0,.1948233))
  BLUE <- "#3B6EA5"; GREY <- "#C8CDD3"; INK <- "#25313C"; ORANGE <- "#D98C3F"
  width <- 12.8; height <- 6.8
  formats <- unique(getOption("ptsd.figure2.formats",c("pdf","png","tiff")))
  if (!length(formats)||!all(formats %in% c("pdf","png","tiff"))) stop("Invalid output formats.")
  parent <- getOption("ptsd.figure2.output_dir",file.path(root,"figures_current"))
  if (!dir.exists(parent) && !dir.create(parent,recursive=TRUE,showWarnings=FALSE)) stop("Cannot create output parent.")
  out <- tempfile(paste0("Figure2_",format(Sys.time(),"%Y%m%d_%H%M%S"),"_"),tmpdir=parent)
  if (!dir.create(out)) stop("Cannot create output folder.")

  draw <- function() {
    graphics::par(omi=c(.35,.22,.25,.28),bg="white",family="sans",col=INK,col.axis=INK,col.lab=INK,
      cex.axis=.83,cex.lab=.90,las=1,tcl=-.25,mgp=c(2.5,.6,0))
    # A has a separate residual note; C has a histogram and aligned candidate strip.
    graphics::layout(matrix(c(1,3,4, 2,3,5),nrow=2,byrow=TRUE),
      widths=c(.92,1.45,1.35),heights=c(3.6,2.1))
    graphics::par(cex=1)

    # A: standardized loadings; residuals are text, never negative bars.
    graphics::par(mar=c(2.6,3.4,3.4,.6),xaxs="i",yaxs="i")
    graphics::plot.new();graphics::plot.window(xlim=c(.4,3.6),ylim=c(0,1.18))
    graphics::rect((1:3)-.29,0,(1:3)+.29,model$standardized_loading,col=BLUE,border=NA)
    graphics::axis(1,at=1:3,labels=model$trait,tick=FALSE,cex.axis=.91,font=2)
    graphics::axis(2,at=c(0,.5,1),labels=c("0","0.5","1.0"))
    graphics::box(bty="l",lwd=.75)
    graphics::text(1:3,model$standardized_loading+.045,sprintf("%.3f",model$standardized_loading),cex=.91)
    graphics::mtext("Standardized loading",side=2,line=2.4,las=0,cex=.90)
    graphics::mtext("A  Common internalizing factor",side=3,line=2.0,adj=0,cex=.91,font=2)
    graphics::mtext("Recorded primary trait model",side=3,line=.65,adj=0,cex=.67,col="#64707C")

    graphics::par(mar=c(1.4,.7,.25,.35),xaxs="i",yaxs="i")
    graphics::plot.new();graphics::plot.window(xlim=c(0,1),ylim=c(0,1))
    graphics::rect(0,.03,1,.98,col="#F5F7F9",border="#DDE2E7",lwd=.65)
    graphics::text(.055,.86,"Standardized residual variance",adj=0,cex=.75,font=2)
    graphics::text(.055,.64,"PTSD: 0.1423",adj=0,cex=.77)
    graphics::text(.055,.47,expression(italic(P)==5.39%*%10^-4),adj=0,cex=.73)
    graphics::text(.055,.28,"MDD: fixed to 0",adj=0,cex=.77)
    graphics::text(.055,.10,"Anxiety: 0.1948",adj=0,cex=.77)

    # B: volcano. Candidate fills follow actual MR-FDR status (CD101 is gray).
    graphics::par(mar=c(4.5,4.0,3.4,.65),xaxs="i",yaxs="i")
    graphics::plot.new();graphics::plot.window(xlim=c(-.65,.65),ylim=c(0,15.2))
    graphics::abline(v=0,col="#E5E8EB",lwd=.65)
    ns<-!d$MR_FDR_significant;sig<-d$MR_FDR_significant
    graphics::points(d$mr_beta[ns],d$MR_minus_log10_P[ns],pch=16,cex=.55,col=grDevices::adjustcolor(GREY,.65))
    graphics::points(d$mr_beta[sig],d$MR_minus_log10_P[sig],pch=16,cex=.70,col=BLUE)
    graphics::points(cd$mr_beta,cd$MR_minus_log10_P,pch=21,cex=1.12,lwd=.85,col=INK,
      bg=ifelse(cd$MR_FDR_significant,BLUE,GREY))
    labels<-data.frame(protein=c("AKT3","CD40","CGREF1","FES","FURIN","KHK","SNX18","UBE2L6","CD101"),
      x=c(-.49,-.25,-.21,.42,-.11,.12,-.39,.22,.18),y=c(5.95,12.45,7.45,14.10,5.0,7.15,4.12,10.3,2.25))
    j<-match(labels$protein,cd$protein)
    graphics::segments(cd$mr_beta[j],cd$MR_minus_log10_P[j],labels$x,labels$y,col="#87939D",lwd=.6)
    graphics::text(labels$x,labels$y,labels$protein,cex=.78,font=2,col=INK)
    graphics::axis(1,at=c(-.5,0,.5),labels=c("-0.5","0","0.5"))
    graphics::axis(2,at=c(0,3,6,9,12,15));graphics::box(bty="l",lwd=.75)
    graphics::mtext("MR effect on common factor (Wald ratio)",side=1,line=2.45,cex=.80)
    graphics::mtext(expression(-log[10](italic(P))),side=2,line=2.5,las=0,cex=.91)
    graphics::mtext("B  Protein-to-common-factor MR",side=3,line=2.0,adj=0,cex=.94,font=2)
    graphics::mtext("1,162 proteins; 51 BH-FDR < 0.05",side=3,line=.65,adj=0,cex=.73,col="#64707C")
    graphics::legend("topright",legend=c("BH-FDR < 0.05","BH-FDR >= 0.05","PTSD-prioritized"),
      pch=c(16,16,21),pt.bg=c(BLUE,GREY,"white"),col=c(BLUE,GREY,INK),pt.cex=c(.75,.75,.95),
      cex=.66,bty="n",inset=c(.015,.14),y.intersp=1.12,x.intersp=.7)

    # C: all protein tests, followed by an aligned strip labeling all 9 candidates.
    xmax<-4.85
    breaks<-seq(0,xmax,by=.15)
    if(tail(breaks,1)<xmax) breaks<-c(breaks,xmax)
    h<-graphics::hist(d$Q_minus_log10_P,breaks=breaks,plot=FALSE,right=FALSE)
    ymax<-max(h$counts)*1.25
    graphics::par(mar=c(.45,4.0,3.4,.65),xaxs="i",yaxs="i")
    graphics::plot.new();graphics::plot.window(xlim=c(0,xmax),ylim=c(0,ymax))
    graphics::rect(head(h$breaks,-1),0,tail(h$breaks,-1),h$counts,col="#9FB8D4",border="white",lwd=.45)
    graphics::axis(2,at=pretty(c(0,max(h$counts)),n=4));graphics::box(bty="l",lwd=.75)
    graphics::abline(v=nominal,lty=2,col="#697783",lwd=.95)
    graphics::abline(v=bonferroni,lty=3,col="#36434E",lwd=.95)
    graphics::text(nominal+.08,ymax*.94,"Nominal\nP = 0.05",adj=c(0,1),cex=.65,col="#596671")
    graphics::text(bonferroni-.08,ymax*.94,"Bonferroni\nP = 0.05/1,162",adj=c(1,1),cex=.65,col="#36434E")
    graphics::mtext("Number of proteins",side=2,line=2.5,las=0,cex=.87)
    graphics::mtext(expression(bold("C  Instrument-level Q"[SNP]*" deviation")),
      side=3,line=2.0,adj=0,cex=.91)
    graphics::mtext("87 nominal; 0 BH-FDR significant",side=3,line=.65,adj=0,cex=.71,col="#64707C")

    graphics::par(mar=c(3.8,4.0,.3,.65),xaxs="i",yaxs="i")
    graphics::plot.new();graphics::plot.window(xlim=c(0,xmax),ylim=c(.5,9.5))
    qc<-cd[order(cd$Q_minus_log10_P,decreasing=TRUE),,drop=FALSE]
    yy<-9:1
    graphics::segments(0,yy,xmax,yy,col="#EEF1F4",lwd=.55)
    graphics::abline(v=nominal,lty=2,col="#697783",lwd=.8)
    graphics::abline(v=bonferroni,lty=3,col="#36434E",lwd=.8)
    graphics::points(qc$Q_minus_log10_P,yy,pch=25,bg=ORANGE,col=INK,cex=.88,lwd=.65)
    graphics::axis(2,at=yy,labels=qc$protein,tick=FALSE,line=-.15,cex.axis=.67)
    graphics::axis(1,at=0:4);graphics::box(bty="l",lwd=.75)
    graphics::mtext(expression(-log[10](italic(P)[Q[SNP]])),side=1,line=2.2,cex=.91)
  }

  export <- function(format) {
    path<-file.path(out,paste0("Figure2_SEM_expanded.",format))
    if(format=="pdf") {
      grDevices::pdf(path,width=width,height=height,family="Helvetica",pointsize=13,
        paper="special",useDingbats=FALSE,bg="white")
    } else if(format=="png") {
      grDevices::png(path,width=width,height=height,units="in",res=600,pointsize=13,bg="white")
    } else {
      grDevices::tiff(path,width=width,height=height,units="in",res=600,pointsize=13,bg="white",compression="lzw")
    }
    on.exit(grDevices::dev.off(),add=TRUE)
    draw();path
  }
  paths<-vapply(formats,export,character(1))
  utils::write.csv(d,file.path(out,"Figure2_plot_data.csv"),row.names=FALSE)
  utils::write.csv(model,file.path(out,"Figure2_trait_model.csv"),row.names=FALSE)
  summary<-data.frame(metric=c(names(counts),"Q_Bonferroni_minus_log10_P"),value=c(unname(counts),bonferroni))
  utils::write.csv(summary,file.path(out,"Figure2_counts.csv"),row.names=FALSE)
  caption<-paste(
    "Figure 2. Common internalizing factor and multivariate protein associations.",
    "(A) Standardized loadings of PTSD, MDD, and anxiety on the recorded primary common factor.",
    "Standardized residual genetic variances are shown as text: PTSD 0.1423 (P = 5.39 x 10^-4), MDD fixed to zero, and anxiety 0.1948.",
    "(B) Protein-to-factor Wald-ratio MR for 1,162 proteins representing 1,152 unique instrument SNPs.",
    "Blue points denote 51 associations with BH-FDR < 0.05; gray points denote the remaining proteins.",
    "Outlined circles and labels identify all nine testable PTSD-prioritized proteins; their fills retain the same FDR color coding.",
    "Eight of the nine candidates passed factor-MR FDR correction; CD101 did not.",
    "(C) Distribution of instrument-level Q_SNP P values across the 1,162 protein tests, with an aligned orange-triangle strip for the nine candidates.",
    "Dashed and dotted lines mark nominal P = 0.05 and Bonferroni P = 0.05/1,162, respectively; the latter is a visual reference, not the primary multiplicity criterion.",
    "Eighty-seven protein tests were nominally significant and none survived BH-FDR correction; proteins sharing a SNP share its Q_SNP statistic.",
    "SIRPA was untestable because its sentinel SNP was absent from MDD summary statistics.",
    "Panels show the recorded primary specification; source-NEFF sensitivity results are reported separately.",
    "Factor-MR results are not independent replication of the PTSD screen, and absence of corrected Q_SNP evidence does not establish equal effects across disorders."
  )
  writeLines(caption,file.path(out,"Figure2_caption.txt"),useBytes=TRUE)
  writeLines(c("Plot-only export of reviewed expanded primary SEM results; no analysis rerun.",
    paste("Input:",normalizePath(input,winslash="/")),paste("Input MD5:",unname(tools::md5sum(input))),
    "Model: constrained primary model, MDD residual fixed at zero; primary loadings unchanged by instrument expansion.",
    "Nine tested candidates; eight factor-MR FDR significant. CD101 is gray in panel B.",
    "SIRPA is absent because MDD data are unavailable; it is not displayed as a null result."),file.path(out,"Figure2_source_note.txt"))
  print(summary,row.names=FALSE)
  cat("\nFigure 2 written to:\n",normalizePath(out,winslash="/"),"\n",sep="")
  cat("PDF + 600-dpi PNG/TIFF are generated by default. No GWAS or SEM rerun.\n")
  invisible(paths)
})()
