# Figure 3 only. Reads complete frozen manuscript tables S3a/S3b/S3c.
# No raw-data read, MR estimation, result reclassification or BH recomputation.

ptsd_figure3_data <- function(snapshot_dir) {
  read <- function(sheet,outcome,expected_n) {
    p <- file.path(snapshot_dir,"current_manuscript_components",paste0(sheet,".csv"))
    if(!file.exists(p)) stop("Missing frozen MR table: ",p)
    d <- read.csv(p,check.names=FALSE,stringsAsFactors=FALSE)
    need <- c("protein","SNP","mr_beta","mr_se","mr_pval","OR_lower","OR_upper","fdr")
    if(!all(need %in% names(d)) || nrow(d)!=expected_n || anyDuplicated(d$protein))
      stop("Unexpected table schema/count: ",sheet)
    numeric <- c("mr_beta","mr_se","mr_pval","OR_lower","OR_upper","fdr")
    if(any(!is.finite(as.matrix(d[numeric]))) || any(d$mr_se<=0) ||
       any(d$mr_pval<=0|d$mr_pval>1|d$fdr<0|d$fdr>1) ||
       any(d$OR_lower<=0|d$OR_upper<=0|d$OR_lower>d$OR_upper))
      stop("Invalid numerical field in ",sheet)
    d$outcome <- outcome
    d
  }
  all <- rbind(read("TableS3a","PTSD",1258L),read("TableS3b","MDD",1330L),
               read("TableS3c","Anxiety",1385L))
  pd <- all[all$outcome=="PTSD",,drop=FALSE]
  pd$nlp <- -log10(pd$mr_pval)
  pd$grp <- ifelse(pd$fdr<0.05,"FDR<0.05","NS")
  display_order <- c("SIRPA","CD40","AKT3","UBE2L6","SNX18",
                     "CGREF1","KHK","FURIN","FES","CD101")
  hits <- pd[pd$fdr<0.05,,drop=FALSE]
  if(!setequal(hits$protein,display_order)) stop("PTSD-prioritized set differs from the reviewed ten proteins")
  b <- all[all$protein %in% display_order,,drop=FALSE]
  b$disorder <- factor(b$outcome,levels=c("PTSD","MDD","Anxiety"))
  b$lo <- log(b$OR_lower);b$hi <- log(b$OR_upper)
  b$sigflag <- b$fdr<0.05
  b$protein <- factor(b$protein,levels=display_order)
  list(volcano=pd,hits=hits,bonferroni=0.05/nrow(pd),
    fdr_line=-log10(max(hits$mr_pval)),forest=b)
}

if(!isTRUE(getOption("ptsd.figure3.load_only",FALSE))) {
  P <- getOption("ptsd.current.paths")
  if(is.null(P)) stop("Source 00_START.R first")
  pkgs <- c("ggplot2","patchwork","ggrepel")
  missing <- pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing)) stop("Missing installed plotting packages: ",paste(missing,collapse=", "))
  suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(ggrepel)})
  d <- ptsd_figure3_data(P$snapshot_dir)
  pd <- d$volcano;hits <- d$hits;b <- d$forest
  bonf_hits <- pd[pd$mr_pval<d$bonferroni,,drop=FALSE]
  blue <- "#3b6ea5";grey <- "#bdbdbd"
  pA <- ggplot(pd,aes(mr_beta,nlp)) +
    geom_point(aes(color=grp),size=1.3,alpha=0.6) +
    scale_color_manual(values=c("FDR<0.05"=blue,"NS"=grey),name=NULL) +
    geom_point(data=bonf_hits,shape=21,fill=NA,color="black",size=3.2,stroke=1) +
    geom_hline(yintercept=d$fdr_line,linetype="dashed",color=blue,linewidth=0.4) +
    geom_hline(yintercept=-log10(d$bonferroni),linetype="dotted",color="black",linewidth=0.4) +
    geom_text_repel(data=hits,aes(label=protein),size=3.4,fontface="bold",seed=42,
      min.segment.length=0,box.padding=0.4) +
    annotate("text",x=max(pd$mr_beta),y=d$fdr_line,label="FDR 0.05",hjust=1,vjust=-0.4,size=3,color=blue) +
    annotate("text",x=max(pd$mr_beta),y=-log10(d$bonferroni),label="Bonferroni",hjust=1,vjust=-0.4,size=3) +
    labs(title="PTSD proteome-wide MR (1,258 proteins)",
      x="MR effect on PTSD\n(log OR per source-scale unit of protein)",y=expression(-log[10](italic(P)))) +
    theme_classic(base_size=12) +
    theme(plot.title=element_text(face="bold",size=13,hjust=0.5),legend.position=c(0.12,0.9))
  dodge <- position_dodge(width=0.6)
  colors <- c(PTSD="#3b6ea5",MDD="#e08214",Anxiety="#4a9b6e")
  pB <- ggplot(b,aes(mr_beta,protein,color=disorder,group=disorder)) +
    geom_vline(xintercept=0,linetype="dashed",color="grey60",linewidth=0.4) +
    geom_errorbarh(aes(xmin=lo,xmax=hi),height=0,position=dodge,linewidth=0.6) +
    geom_point(aes(fill=disorder),position=dodge,size=2.4,shape=21,color="white",stroke=0.3) +
    # Keep all rows in this layer so significant-point outlines are not dodged
    # to different positions when another disorder is non-significant/missing.
    geom_point(aes(alpha=as.numeric(sigflag)),position=dodge,size=2.6,shape=21,
      fill=NA,color="black",stroke=0.9) +
    scale_alpha_identity(guide="none") +
    scale_color_manual(values=colors,name=NULL) +
    scale_fill_manual(values=colors,guide="none") +
    labs(title="Per-disorder MR estimates for prioritized proteins",
      x="MR effect\n(log OR per source-scale unit of protein)",y=NULL,
      caption="Dark-outlined point = FDR<0.05 in that disorder") +
    theme_classic(base_size=12) +
    theme(plot.title=element_text(face="bold",size=11.5,hjust=0.5),legend.position="bottom")
  fig <- (pA|pB)+plot_layout(widths=c(1,1))+plot_annotation(tag_levels="A") &
    theme(plot.tag=element_text(face="bold",size=14))
  ggsave(file.path(P$report_dir,"Figure3_current.pdf"),fig,width=14,height=6)
  ggsave(file.path(P$report_dir,"Figure3_current.png"),fig,width=14,height=6,dpi=300)
  write.csv(b,file.path(P$report_dir,"Figure3_forest_plot_values.csv"),row.names=FALSE)
  cat("Figure 3 written from frozen complete MR tables; source-scale labels, no association rerun.\n")
}
