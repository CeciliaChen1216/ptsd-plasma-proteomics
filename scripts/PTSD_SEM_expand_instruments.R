# Expand the Genomic SEM instrument universe with an auditable reference/QC step.
# Uses configured raw GWAS, the 1000G EUR panel and reviewed LDSC caches.
# No packages are installed and no old analysis files are overwritten.
# Run after 00_START.R: ptsd_run("sem_expand")
# Optional: options(ptsd.sem.expand.audit_only=TRUE) for preparation/audit only.
# Optional: options(ptsd.sem.expand.resume_dir="existing extraction directory").
# Resume requires an explicit option; the release defaults to fresh extraction.
# A new timestamped result directory is created under the configured project.
# Genome-wide LDSC is NOT restricted to the target instrument SNPs or rerun here.
# The SNP step covers ALL instruments in the existing pQTL selection, not only
# significant PTSD candidates. Palindromes require independently checked EAFs.
# Compatibility fix 2026-10-18: reference coordinates are optional; missing
# chromosome/position metadata are filled from the audited GRCh37 EUR .bim.
# Compatibility fix 2026-10-19: recognize PGC FRQ_U_<N>/FRQ_A_<N> EAFs;
# retain original nonpalindromic reference SNPs using concordant GRCh37 GWAS
# coordinates when absent from EUR. A reviewed failed extraction can be resumed.

semx_col <- function(nms, aliases, required=TRUE) {
  key <- toupper(sub("^#+", "", trimws(sub("^\ufeff", "", nms))))
  hit <- match(toupper(aliases), key, nomatch=0L)
  hit <- hit[hit>0L]
  if (!length(hit)) {
    if (required) stop("Missing column; expected one of: ",paste(aliases,collapse=", "),
                       ". Available columns: ",paste(nms,collapse=", "))
    return(NA_character_)
  }
  nms[hit[1L]]
}
semx_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

semx_read_reference <- function(path,ids) {
  h<-names(data.table::fread(path,nrows=0L,showProgress=FALSE))
  required<-list(SNP=c("SNP","RSID","RS_ID","VARIANT_ID","ID"),MAF="MAF",A1="A1",A2="A2")
  map<-vapply(required,function(a) semx_col(h,a),character(1))
  map<-c(map,CHR=semx_col(h,c("CHR","CHROM","CHROMOSOME"),FALSE),
             BP=semx_col(h,c("BP","POS","POSITION","BASE_PAIR_LOCATION"),FALSE))
  present<-!is.na(map)
  rr<-data.table::fread(path,select=unname(map[present]),data.table=FALSE,showProgress=FALSE)
  names(rr)<-names(map)[present]
  total<-nrow(rr)
  rr<-rr[rr$SNP%in%ids,,drop=FALSE]
  for(nm in names(map)[!present]) rr[[nm]]<-rep(NA_real_,nrow(rr))
  rr<-rr[,c("SNP","CHR","BP","MAF","A1","A2"),drop=FALSE]
  rr$A1<-toupper(rr$A1);rr$A2<-toupper(rr$A2)
  rr$CHR<-semx_num(sub("^chr","",rr$CHR,ignore.case=TRUE))
  rr$BP<-semx_num(rr$BP);rr$MAF<-semx_num(rr$MAF)
  list(data=rr,total_rows=total,header=h,mapping=map)
}

semx_reference_coordinates <- function(ref,bim) {
  ix<-match(ref$SNP,bim$SNP)
  audit<-data.frame(SNP=ref$SNP,original_CHR=ref$CHR,original_BP=ref$BP,
    CHR_source=ifelse(is.na(ref$CHR),"unavailable","original_reference"),
    BP_source=ifelse(is.na(ref$BP),"unavailable","original_reference"),stringsAsFactors=FALSE)
  chr_fill<-is.na(ref$CHR)&!is.na(ix)&is.finite(bim$CHR[ix])
  bp_fill<-is.na(ref$BP)&!is.na(ix)&is.finite(bim$BP[ix])
  ref$CHR[chr_fill]<-bim$CHR[ix[chr_fill]]
  ref$BP[bp_fill]<-bim$BP[ix[bp_fill]]
  audit$CHR_source[chr_fill]<-"audited_1000G_EUR_GRCh37"
  audit$BP_source[bp_fill]<-"audited_1000G_EUR_GRCh37"
  audit$resolved_CHR<-ref$CHR;audit$resolved_BP<-ref$BP
  list(data=ref,audit=audit)
}

semx_frequency_columns <- function(nms,trait) {
  keys<-toupper(nms)
  aliases<-c("EAF","EFFECT_ALLELE_FREQUENCY","A1FREQ","A1_FREQ","FRQ","FREQ")
  take<-match(aliases,keys,nomatch=0L);take<-take[take>0L]
  # In PGC daner files FRQ_A is A1 frequency in cases and FRQ_U in controls;
  # the optional numeric suffix describes the contributing sample count.
  controls<-grep("^FRQ_U(_[0-9]+)?$",keys)
  cases<-grep("^FRQ_A(_[0-9]+)?$",keys)
  if(length(controls)>1L||length(cases)>1L)
    stop("Multiple case/control EAF columns: inspect ",paste(nms,collapse=", "))
  take<-unique(c(take,controls,cases,if(trait=="PTSD") which(keys=="MAF")))
  nms[take]
}

semx_gwas_coordinates <- function(ref,raws) {
  audited<-data.frame(SNP=ref$SNP,previous_CHR=ref$CHR,previous_BP=ref$BP,
    coordinate_action="reference_coordinates_present",stringsAsFactors=FALSE)
  need<-which(is.na(ref$CHR)|is.na(ref$BP))
  for(trait in c("MDD","ANX")) {
    x<-raws[[trait]];sp<-semx_spec(x,trait)
    id<-as.character(x[[sp$SNP]]);dup<-id[duplicated(id)|duplicated(id,fromLast=TRUE)]
    ix<-match(ref$SNP,id)
    chr<-if(!is.na(sp$CHR)) semx_num(sub("^chr","",x[[sp$CHR]][ix],ignore.case=TRUE)) else rep(NA_real_,nrow(ref))
    bp<-if(!is.na(sp$BP)) semx_num(x[[sp$BP]][ix]) else rep(NA_real_,nrow(ref))
    ea<-toupper(x[[sp$A1]][ix]);oa<-toupper(x[[sp$A2]][ix])
    pair<-((ea==ref$A1&oa==ref$A2)|(ea==ref$A2&oa==ref$A1)|
      (semx_comp(ea)==ref$A1&semx_comp(oa)==ref$A2)|
      (semx_comp(ea)==ref$A2&semx_comp(oa)==ref$A1))&!ref$SNP%in%dup
    pair[is.na(pair)]<-FALSE
    audited[[paste0(trait,"_CHR")]]<-chr
    audited[[paste0(trait,"_BP")]]<-bp
    audited[[paste0(trait,"_allele_pair_matches")]]<-pair
  }
  for(i in need) {
    ch<-c(audited$MDD_CHR[i],audited$ANX_CHR[i]);bp<-c(audited$MDD_BP[i],audited$ANX_BP[i])
    good<-all(is.finite(c(ch,bp)))&&all(ch%in%1:22)&&all(bp>0)&&all(bp==floor(bp))&&
      ch[1]==ch[2]&&bp[1]==bp[2]&&audited$MDD_allele_pair_matches[i]&&audited$ANX_allele_pair_matches[i]&&
      (is.na(ref$CHR[i])||ref$CHR[i]==ch[1])&&(is.na(ref$BP[i])||ref$BP[i]==bp[1])
    if(good) {
      ref$CHR[i]<-ch[1];ref$BP[i]<-bp[1]
      audited$coordinate_action[i]<-"filled_from_concordant_MDD_ANX_GRCh37"
    } else audited$coordinate_action[i]<-"coordinates_unresolved_or_discordant"
  }
  audited$resolved_CHR<-ref$CHR;audited$resolved_BP<-ref$BP
  list(data=ref,audit=audited)
}

semx_resume <- function(root,required,q) {
  path<-getOption("ptsd.sem.expand.resume_dir","")
  if(!nzchar(path)) return(NULL)
  path<-normalizePath(path,winslash="/",mustWork=TRUE)
  needed<-c("input_manifest.csv","reference_target_audit.csv","pQTL_harmonization_audit.csv",
    paste0(c("PTSD","MDD","ANX"),"_raw_instrument_rows.csv"),"EUR_instrument_frequencies.afreq")
  if(any(!file.exists(file.path(path,needed)))) stop("Resume directory lacks required extraction files: ",path)
  m<-data.table::fread(file.path(path,"input_manifest.csv"),data.table=FALSE)
  norm<-function(x) tolower(gsub("\\\\","/",x))
  ix<-match(norm(required),norm(m$path));now<-file.info(required)
  if(anyNA(ix)||any(now$size!=m$bytes[ix])) stop("Resume inputs no longer match the recorded file paths/sizes.")
  recorded<-as.POSIXct(m$modified[ix],format="%Y-%m-%dT%H:%M:%OSZ",tz="UTC")
  if(anyNA(recorded)||any(abs(as.numeric(now$mtime)-as.numeric(recorded))>2))
    stop("Resume input modification times differ; inspect before reusing target rows.")
  oldq<-data.table::fread(file.path(path,"pQTL_harmonization_audit.csv"),data.table=FALSE)
  ix<-match(q$protein,oldq$protein)
  if(nrow(q)!=nrow(oldq)||anyNA(ix)||any(q$SNP!=oldq$SNP[ix])||
     any(abs(q$beta-oldq$beta[ix])>1e-12)||any(abs(q$se-oldq$se[ix])>1e-12))
    stop("Resume pQTL instrument selection differs from the recorded selection.")
  list(dir=path,files=needed)
}

semx_comp <- function(x) chartr("ACGT","TGCA",toupper(x))
semx_pal <- function(a,b) paste0(a,b) %in% c("AT","TA","CG","GC")
semx_align <- function(ea,oa,eaf,r1,r2,ref_eaf,tol=.10,margin=.10) {
  ea<-toupper(ea);oa<-toupper(oa);r1<-toupper(r1);r2<-toupper(r2)
  fail<-function(reason) list(ok=FALSE,sign=NA_real_,action=reason,aligned_eaf=NA_real_)
  if(anyNA(c(ea,oa,r1,r2)) || !all(c(ea,oa,r1,r2)%in%c("A","C","G","T")) ||
     ea==oa || r1==r2) return(fail("invalid_or_non_SNV_alleles"))
  if(semx_pal(ea,oa)) {
    if(!semx_pal(r1,r2) || !setequal(c(ea,oa),c(r1,r2))) return(fail("allele_pair_mismatch"))
    if(!is.finite(eaf)||eaf<=0||eaf>=1||!is.finite(ref_eaf)||ref_eaf<=0||ref_eaf>=1)
      return(fail("palindrome_missing_effect_allele_frequency"))
    if(min(ref_eaf,1-ref_eaf)>.42) return(fail("palindrome_reference_MAF_above_0.42"))
    s0<-if(ea==r1) 1 else -1
    f0<-if(s0==1) eaf else 1-eaf
    dist<-c(abs(f0-ref_eaf),abs(1-f0-ref_eaf)); win<-which.min(dist)
    if(dist[win]>tol || dist[3L-win]<=tol || diff(sort(dist))<margin)
      return(fail("palindrome_frequency_ambiguous_or_discordant"))
    sg<-if(win==1L) s0 else -s0
    return(list(ok=TRUE,sign=sg,action=if(win==1L) "palindrome_frequency_resolved" else
                  "palindrome_frequency_resolved_strand_flip",
                aligned_eaf=if(sg==1) eaf else 1-eaf))
  }
  sg<-if(ea==r1&&oa==r2) 1 else if(ea==r2&&oa==r1) -1 else
      if(semx_comp(ea)==r1&&semx_comp(oa)==r2) 1 else
      if(semx_comp(ea)==r2&&semx_comp(oa)==r1) -1 else NA_real_
  if(!is.finite(sg)) return(fail("allele_pair_mismatch"))
  list(ok=TRUE,sign=sg,action=if(ea%in%c(r1,r2)) "allele_pair_aligned" else "strand_complement_aligned",
       aligned_eaf=if(is.finite(eaf)) if(sg==1) eaf else 1-eaf else NA_real_)
}

semx_read_targets <- function(path, ids, label, chunk=100000L) {
  con<-if(grepl("[.]gz$",path,ignore.case=TRUE)) gzfile(path,"rt") else file(path,"rt")
  on.exit(close(con),add=TRUE)
  repeat {
    header<-readLines(con,n=1L,warn=FALSE)
    if(!length(header)) stop("No header in ",path)
    if(nzchar(trimws(header))&&!startsWith(header,"##")) break
  }
  nms<-strsplit(trimws(header),"[[:space:]]+")[[1L]]
  id_col<-semx_col(nms,c("SNP","RSID","RS_ID","VARIANT_ID","ID","MARKERNAME"))
  at<-match(id_col,nms)
  pat<-paste0("^[[:space:]]*(?:\\S+\\s+){",at-1L,"}(\\S+).*$")
  sep<-if(grepl("\t",header,fixed=TRUE)) "\t" else "auto"
  found<-list();total<-0L;k<-0L
  repeat {
    lines<-readLines(con,n=chunk,warn=FALSE)
    if(!length(lines)) break
    total<-total+length(lines);k<-k+1L
    snps<-sub(pat,"\\1",lines,perl=TRUE)
    hit<-data.table::chmatch(snps,ids,nomatch=0L)>0L
    if(any(hit)) found[[length(found)+1L]]<-data.table::fread(
      text=paste(c(header,lines[hit]),collapse="\n"),sep=sep,header=TRUE,
      colClasses="character",data.table=FALSE,check.names=FALSE,showProgress=FALSE)
    if(k%%10L==0L) {cat(label,": scanned ",format(total,big.mark=",")," rows\n",sep="");flush.console()}
  }
  cat(label,": ",format(total,big.mark=",")," rows read; target rows=",
      sum(vapply(found,nrow,integer(1))),"\n",sep="")
  if(!length(found)) {
    out<-as.data.frame(setNames(replicate(length(nms),character(),simplify=FALSE),nms),check.names=FALSE)
  } else out<-as.data.frame(data.table::rbindlist(found,use.names=TRUE,fill=TRUE))
  attr(out,"rows_scanned")<-total
  out
}

semx_spec <- function(x,trait) {
  n<-names(x)
  eff<-semx_col(n,c("BETA","B","LOG_OR","LOGOR","LOG_ODDS"),FALSE)
  kind<-"beta"
  if(is.na(eff)) {eff<-semx_col(n,c("OR","ODDS_RATIO"));kind<-"OR"}
  eafs<-semx_frequency_columns(n,trait)
  # PTSD MAF is the known effect-allele frequency renamed from .ma$freq, as in
  # the reviewed preprocessing. Generic MAF in other files does not identify EA.
  list(SNP=semx_col(n,c("SNP","RSID","RS_ID","VARIANT_ID","ID","MARKERNAME")),
       A1=semx_col(n,c("A1","EA","EFFECT_ALLELE","ALLELE1")),
       A2=semx_col(n,c("A2","OA","OTHER_ALLELE","ALLELE0","NEA")),
       effect=eff,kind=kind,SE=semx_col(n,c("SE","STANDARD_ERROR","STDERR","SEBETA","SE_BETA")),
       P=semx_col(n,c("P","PVAL","P_VALUE","PVALUE")),
       EAF=if(length(eafs)) eafs[1L] else NA_character_,INFO=semx_col(n,c("INFO","INFO_SCORE","IMPUTATION_INFO"),FALSE),
       CHR=semx_col(n,c("CHR","CHROM","CHROMOSOME"),FALSE),
       BP=semx_col(n,c("BP","POS","BASE_PAIR_LOCATION"),FALSE))
}

semx_prepare_trait <- function(x,trait,ref,ids) {
  sp<-semx_spec(x,trait)
  sid<-as.character(x[[sp$SNP]])
  duplicates<-unique(sid[duplicated(sid)|duplicated(sid,fromLast=TRUE)])
  audit<-data.frame(SNP=ids,trait=trait,raw_present=ids%in%sid,
    reason="raw_SNP_absent",raw_A1=NA_character_,raw_A2=NA_character_,raw_beta=NA_real_,
    raw_SE=NA_real_,raw_P=NA_character_,raw_EAF=NA_real_,EAF_column=sp$EAF,
    INFO=NA_real_,alignment=NA_character_,sign=NA_real_,stringsAsFactors=FALSE)
  keep<-list()
  freq_cols<-semx_frequency_columns(names(x),trait)
  for(j in seq_along(ids)) {
    id<-ids[j];i<-match(id,sid);r<-match(id,ref$SNP)
    if(is.na(i)) next
    audit$reason[j]<-"duplicate_raw_SNP_ID"
    if(id%in%duplicates) next
    ea<-toupper(x[[sp$A1]][i]);oa<-toupper(x[[sp$A2]][i])
    e<-semx_num(x[[sp$effect]][i]);b<-if(sp$kind=="OR") if(is.finite(e)&&e>0) log(e) else NA_real_ else e
    se<-semx_num(x[[sp$SE]][i]);pv<-as.character(x[[sp$P]][i]);pn<-semx_num(pv)
    ef<-if(!is.na(sp$EAF)) semx_num(x[[sp$EAF]][i]) else NA_real_
    inf<-if(!is.na(sp$INFO)) semx_num(x[[sp$INFO]][i]) else NA_real_
    audit[j,c("raw_A1","raw_A2","raw_beta","raw_SE","raw_P","raw_EAF","INFO")]<-
      list(ea,oa,b,se,pv,ef,inf)
    audit$reason[j]<-"reference_unavailable_or_failed_QC"
    if(is.na(r)) next
    audit$reason[j]<-"invalid_or_zero_effect_SE_or_P"
    if(!is.finite(b)||b==0||!is.finite(se)||se<=0||!is.finite(pn)||pn<=0||pn>1) next
    audit$reason[j]<-"INFO_below_0.6_or_missing_in_existing_column"
    if(!is.na(sp$INFO)&&(!is.finite(inf)||inf<.6)) next
    audit$reason[j]<-"GRCh37_coordinate_mismatch"
    if(!is.na(sp$BP)) {
      bp<-semx_num(x[[sp$BP]][i]); if(is.finite(bp)&&bp!=ref$BP[r]) next
    }
    if(!is.na(sp$CHR)) {
      ch<-semx_num(sub("^chr","",x[[sp$CHR]][i],ignore.case=TRUE))
      if(is.finite(ch)&&ch!=ref$CHR[r]) next
    }
    a<-semx_align(ea,oa,ef,ref$A1[r],ref$A2[r],ref$ref_EAF[r])
    audit$alignment[j]<-a$action; audit$sign[j]<-a$sign
    audit$reason[j]<-a$action
    if(!a$ok) next
    row<-x[i,,drop=FALSE]
    row[[sp$A1]]<-ref$A1[r];row[[sp$A2]]<-ref$A2[r]
    row[[sp$effect]]<-sprintf("%.17g",if(sp$kind=="OR") exp(b*a$sign) else b*a$sign)
    for(f in freq_cols) {
      q<-semx_num(row[[f]])
      if(is.finite(q)) row[[f]]<-sprintf("%.17g",if(a$sign==1) q else 1-q)
    }
    audit$reason[j]<-"included"
    keep[[length(keep)+1L]]<-row
  }
  out<-if(length(keep)) as.data.frame(data.table::rbindlist(keep,use.names=TRUE)) else x[FALSE,,drop=FALSE]
  # Guard the package's median-based OR/log-beta interpretation on small inputs.
  vals<-semx_num(out[[sp$effect]])
  if(length(vals)&&round(median(vals,na.rm=TRUE))!=if(sp$kind=="OR") 1 else 0)
    stop(trait,": median effect would be misclassified by sumstats; inspect raw data.")
  list(data=out,audit=audit,spec=sp)
}

semx_baseline <- function(root,P) {
  paths<-list.files(file.path(root,"results"),pattern="^LDSCoutput_stage0[.]rds$",recursive=TRUE,full.names=TRUE)
  paths<-paths[!grepl("sem_instrument_expansion_",paths,fixed=TRUE)]
  frozen<-P$sem_baseline_snapshot
  if(is.null(frozen)) stop("Missing pinned baseline snapshot setting; source this release 00_START.R first.")
  if(!file.exists(frozen)) stop("Missing the reviewed historical baseline snapshot.")
  reviewed<-data.table::fread(frozen,data.table=FALSE)
  found<-list();status<-list()
  for(path in paths) {
    od<-dirname(path);f<-file.path(od,"parameters.txt")
    pars<-tryCatch(dget(f),error=function(e) NULL)
    kind<-if(is.list(pars)&&is.numeric(pars$sample.prev)&&length(pars$sample.prev)==3L&&
              isTRUE(all.equal(as.numeric(pars$sample.prev),c(137136/(137136+1085746),.5,.5),tolerance=1e-8))) "primary" else
          if(is.list(pars)&&is.numeric(pars$sample.prev)&&length(pars$sample.prev)==3L&&
              isTRUE(all.equal(as.numeric(pars$sample.prev),c(.5,.5,.5),tolerance=1e-8))) "NEFF" else "other"
    required<-file.path(od,c("instrument_sumstats.csv","userGWAS_raw.csv","stage12_factorMR_results.csv","model_primary.txt"))
    ok<-kind!="other"&&all(file.exists(required))&&
      identical(pars$trait_order,c("PTSD","MDD","ANX"))&&
      isTRUE(all.equal(as.numeric(pars$N),c(NA_real_,829250,390020)))&&
      isTRUE(all.equal(as.numeric(pars$population.prev),c(.07,.15,.10)))
    if(ok) {
      v<-tryCatch(data.table::fread(required[3],data.table=FALSE),error=function(e) NULL)
      if(is.null(v)||!all(c("protein","mr_beta","Q_pval")%in%names(v))) ok<-FALSE
      if(ok&&kind=="primary") {
        ix<-match(reviewed$protein,v$protein)
        ok<-nrow(v)==nrow(reviewed)&&!anyNA(ix)&&
          all(abs(v$mr_beta[ix]-reviewed$mr_beta)<1e-8)&&
          all(abs(v$Q_pval[ix]-reviewed$Q_pval)<1e-7)
      }
    }
    status[[length(status)+1L]]<-data.frame(kind=kind,directory=od,eligible=ok)
    if(ok) found[[length(found)+1L]]<-list(kind=kind,dir=od,pars=pars,results=v)
  }
  if(!length(found)) stop("No reviewed LDSC cache with instrument_sumstats.csv and userGWAS_raw.csv found.")
  choose<-function(kind) {
    zz<-Filter(function(v) v$kind==kind,found)
    explicit<-getOption(paste0("ptsd.sem.expand.",tolower(kind),"_dir"),
      if(kind=="primary") P$sem_primary_cache else P$sem_neff_cache)
    if(nzchar(explicit)) zz<-Filter(function(v) normalizePath(v$dir,winslash="/")==normalizePath(explicit,winslash="/"),zz)
    if(!length(zz)) return(NULL)
    if(length(zz)>1L) {
      hashes<-vapply(zz,function(v) unname(tools::md5sum(file.path(v$dir,"LDSCoutput_stage0.rds"))),character(1))
      if(length(unique(hashes))>1L) stop("Different ",kind," LDSC caches found. Set options(ptsd.sem.expand.",
        tolower(kind),"_dir='the reviewed run directory').")
    }
    zz[[1L]]
  }
  primary<-choose("primary")
  if(is.null(primary)) stop("The recorded primary LDSC cache could not be verified against the reviewed 391-protein results.")
  list(primary=primary,NEFF=choose("NEFF"),inventory=do.call(rbind,status),reviewed=reviewed)
}

semx_same_ss <- function(a,b,ids,cols,tol=1e-9) {
  ia<-match(ids,a$SNP);ib<-match(ids,b$SNP);ok<-!is.na(ia)&!is.na(ib)
  # Reference SNP variance and allele orientation are part of the same input.
  if(all(c("MAF","A1","A2")%in%names(a))&&all(c("MAF","A1","A2")%in%names(b))) {
    cols<-unique(c("MAF",cols))
    ok<-ok&!is.na(a$A1[ia])&!is.na(b$A1[ib])&!is.na(a$A2[ia])&!is.na(b$A2[ib])&
      a$A1[ia]==b$A1[ib]&a$A2[ia]==b$A2[ib]
  }
  for(col in cols) {
    x<-semx_num(a[[col]][ia]);y<-semx_num(b[[col]][ib])
    ok<-ok&is.finite(x)&is.finite(y)&abs(x-y)<=tol*pmax(1,abs(y))
  }
  ok
}

semx_run <- function() {
  P<-getOption("ptsd.current.paths")
  if(is.null(P)) stop("Source 00_START.R first, then use ptsd_run('sem_expand').")
  root<-P$project_root
  dir.create(file.path(root,"results"),recursive=TRUE,showWarnings=FALSE)
  out<-tempfile(paste0("sem_instrument_expansion_",format(Sys.time(),"%Y%m%d_%H%M%S"),"_"),tmpdir=file.path(root,"results"))
  if(!dir.create(out)) stop("Cannot create output directory.")
  wd<-getwd();sink_depth<-sink.number();msg_depth<-sink.number(type="message")
  con<-file(file.path(out,"RUN.log"),"wt");sink(con,split=TRUE);sink(con,type="message")
  on.exit({if(sink.number(type="message")!=msg_depth) sink(type="message")
           while(sink.number()>sink_depth) sink();close(con);setwd(wd)
           cat("\nOutput directory: ",out,"\n",sep="")},add=TRUE)
  tryCatch({
    needed<-c("data.table","GenomicSEM","R.utils")
    missing<-needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
    if(length(missing)) stop("Missing installed package(s): ",paste(missing,collapse=", "))
    if(!"Q_SNP"%in%names(formals(GenomicSEM::userGWAS))) stop("Installed userGWAS does not expose Q_SNP.")
    setwd(out)
    writeLines(capture.output(sessionInfo()),"sessionInfo.txt")
    writeLines(capture.output(packageDescription("GenomicSEM")),"GenomicSEM_description.txt")
    for(fn in c("sumstats","userGWAS")) writeLines(capture.output(dput(getExportedValue("GenomicSEM",fn))),paste0("installed_",fn,".R.txt"))
    files<-c(PTSD=P$ptsd_ma,MDD=P$mdd_gwas,ANX=P$anx_gwas)
    required<-c(files,P$pqtl_rds,P$sumstats_reference,P$plink2,paste0(P$ld_prefix,c(".bed",".bim",".fam")))
    if(any(!file.exists(required))) stop("Missing input: ",paste(required[!file.exists(required)],collapse="; "))
    if(basename(files["MDD"])!="daner_pgc_mdd_no23andMe-noUKBB_eur_hg19_v3.49.24.11.neff.gz" ||
       basename(files["ANX"])!="ANX_2026_daner_fullANX_v12_woUTAH_11022026.gz") stop("Unexpected outcome release filename.")
    data.table::fwrite(data.frame(path=required,bytes=file.info(required)$size,modified=file.info(required)$mtime),"input_manifest.csv")
    cat("STEP 1/6: Resolve reviewed LDSC caches and the full pQTL instrument universe\n")
    base<-semx_baseline(root,P);data.table::fwrite(base$inventory,"LDSC_cache_inventory.csv")
    if(is.null(base$NEFF)) stop("A reviewed source-NEFF LDSC cache is required to rebuild both release specifications. Run sem_neff or configure its cache directory.")
    q<-as.data.frame(readRDS(P$pqtl_rds));must<-c("protein","SNP","beta","se","effect_allele","other_allele")
    if(!all(must%in%names(q))) stop("pQTL selection requires: ",paste(must,collapse=", "))
    if(anyDuplicated(q$protein)) stop("pQTL selection contains multiple instruments per protein.")
    q$effect_allele<-toupper(q$effect_allele);q$other_allele<-toupper(q$other_allele)
    q$beta<-semx_num(q$beta);q$se<-semx_num(q$se)
    q$F_stat<-q$beta^2/q$se^2
    q$valid_pqtl<-!is.na(q$SNP)&nzchar(q$SNP)&is.finite(q$beta)&q$beta!=0&is.finite(q$se)&q$se>0&q$F_stat>10
    ids<-unique(q$SNP[q$valid_pqtl]);ids<-ids[grepl("^rs[0-9]+$",ids)]
    if(!length(ids)) stop("No valid selected instruments.")
    cat("Selected proteins=",nrow(q),"; valid rsIDs=",length(ids),"\n",sep="")

    resume<-semx_resume(root,required,q)
    cat("STEP 2/6: Audit old reference and obtain independent EUR frequencies\n")
    if(!is.null(resume)) {
      cat("RESUME: using saved target rows and EUR frequencies from\n",resume$dir,"\n",sep="")
      writeLines(resume$dir,"RESUMED_FROM.txt")
      ref<-data.table::fread(file.path(resume$dir,"reference_target_audit.csv"),data.table=FALSE)
      if(!all(c("SNP","CHR","BP","MAF","A1","A2","old_reference","ref_EAF","reason")%in%names(ref)))
        stop("Resume reference audit schema differs.")
      old_ids<-ref$SNP[ref$old_reference]
      prior_summary<-readLines(file.path(resume$dir,"reference_summary.txt"),warn=FALSE)
      old_reference_n<-semx_num(sub("^Old reference rows: *","",grep("^Old reference rows:",prior_summary,value=TRUE)))
      for(f in c("reference_column_mapping.txt","reference_coordinate_audit.csv","EUR_instrument_frequencies.afreq",
                 "EUR_instrument_frequencies.log","PLINK_command.txt","PLINK_console.txt")) {
        if(file.exists(file.path(resume$dir,f))) file.copy(file.path(resume$dir,f),file.path(out,f))
      }
    } else {
    legacy_ref<-semx_read_reference(P$sumstats_reference,ids)
    rr<-legacy_ref$data;old_reference_n<-legacy_ref$total_rows
    writeLines(c(paste("Original reference columns:",paste(legacy_ref$header,collapse=", ")),
      capture.output(dput(legacy_ref$mapping))),"reference_column_mapping.txt")
    dup<-rr$SNP[duplicated(rr$SNP)|duplicated(rr$SNP,fromLast=TRUE)]
    old_ids<-unique(rr$SNP)
    rr<-rr[!rr$SNP%in%dup,,drop=FALSE]
    b<-data.table::fread(paste0(P$ld_prefix,".bim"),header=FALSE,select=c(1,2,4,5,6),
                         colClasses="character",data.table=FALSE,showProgress=FALSE)
    names(b)<-c("CHR","SNP","BP","B1","B2");b<-b[b$SNP%in%ids,,drop=FALSE]
    bd<-b$SNP[duplicated(b$SNP)|duplicated(b$SNP,fromLast=TRUE)]
    b<-b[!b$SNP%in%bd,,drop=FALSE];b$CHR<-semx_num(b$CHR);b$BP<-semx_num(b$BP)
    anchors<-c(rs1894401=91429042,rs2071410=91420940)
    ai<-match(names(anchors),b$SNP)
    if(anyNA(ai)||any(b$CHR[ai]!=15)||any(b$BP[ai]!=unname(anchors))) stop("EUR reference anchor positions are not the audited GRCh37 positions.")
    rc<-semx_reference_coordinates(rr,b)
    rr<-rc$data;data.table::fwrite(rc$audit,"reference_coordinate_audit.csv")
    cat("Old reference coordinates filled from EUR: CHR=",sum(rc$audit$CHR_source=="audited_1000G_EUR_GRCh37"),
        "; BP=",sum(rc$audit$BP_source=="audited_1000G_EUR_GRCh37"),"\n",sep="")
    fam<-data.table::fread(paste0(P$ld_prefix,".fam"),header=FALSE,select=c(1,2),data.table=FALSE,showProgress=FALSE)
    if(nrow(fam)!=503L) stop("Expected the reviewed 503-participant EUR reference.")
    writeLines(b$SNP,"EUR_frequency_target_ids.txt")
    args<-c("--bfile",P$ld_prefix,"--extract",file.path(out,"EUR_frequency_target_ids.txt"),"--freq",
            "--threads","2","--out",file.path(out,"EUR_instrument_frequencies"))
    writeLines(paste(shQuote(P$plink2),paste(vapply(args,shQuote,character(1)),collapse=" ")),"PLINK_command.txt")
    console<-suppressWarnings(system2(P$plink2,args=vapply(args,shQuote,character(1)),stdout=TRUE,stderr=TRUE))
    writeLines(console,"PLINK_console.txt");exit_status<-attr(console,"status")
    if(!is.null(exit_status)&&exit_status!=0L) stop("PLINK frequency extraction failed.")
    af<-data.table::fread("EUR_instrument_frequencies.afreq",data.table=FALSE)
    if(!all(c("ID","REF","ALT","ALT_FREQS")%in%names(af))) stop("Unexpected PLINK frequency schema.")
    af$ALT_FREQS<-semx_num(af$ALT_FREQS)
    af<-af[!duplicated(af$ID)&!duplicated(af$ID,fromLast=TRUE)&af$REF%in%c("A","C","G","T")&
      af$ALT%in%c("A","C","G","T")&is.finite(af$ALT_FREQS)&af$ALT_FREQS>0&af$ALT_FREQS<1,,drop=FALSE]
    add<-b[!b$SNP%in%c(old_ids,dup),c("SNP","CHR","BP"),drop=FALSE]
    ix<-match(add$SNP,af$ID);add$MAF<-pmin(af$ALT_FREQS[ix],1-af$ALT_FREQS[ix]);add$A1<-af$ALT[ix];add$A2<-af$REF[ix]
    ref<-rbind(rr,add[,names(rr),drop=FALSE]);ref$old_reference<-ref$SNP%in%old_ids
    fi<-match(ref$SNP,af$ID)
    ref$ref_EAF<-ifelse(ref$A1==af$ALT[fi]&ref$A2==af$REF[fi],af$ALT_FREQS[fi],
      ifelse(ref$A1==af$REF[fi]&ref$A2==af$ALT[fi],1-af$ALT_FREQS[fi],NA_real_))
    bi<-match(ref$SNP,b$SNP)
    ref$reason<-"included"
    mismatch<-!is.na(bi)&(ref$CHR!=b$CHR[bi]|ref$BP!=b$BP[bi]);mismatch[is.na(mismatch)]<-FALSE
    ref$reason[mismatch]<-"reference_coordinate_mismatch"
    }

    cat("STEP 3/6: Read target variants from the three original GWAS and harmonize\n")
    preps<-list();raws<-list()
    for(trait in names(files)) {
      if(!is.null(resume)) {
        x<-data.table::fread(file.path(resume$dir,paste0(trait,"_raw_instrument_rows.csv")),
          colClasses="character",data.table=FALSE,check.names=FALSE)
        cat(trait,": reused ",nrow(x)," previously extracted target rows\n",sep="")
      } else x<-semx_read_targets(files[trait],ids,trait)
      if(trait=="PTSD"&&is.null(resume)) {
        expected<-c("SNP","A1","A2","freq","b","se","p","N")
        if(!all(expected%in%names(x))) stop("Unexpected PTSD .ma schema.")
        names(x)[match(expected,names(x))]<-c("SNP","A1","A2","MAF","BETA","SE","P","N")
      }
      data.table::fwrite(x,paste0(trait,"_raw_instrument_rows.csv"));raws[[trait]]<-x
    }
    # A nonpalindromic SNP does not require LD for single-SNP SEM. Missing
    # coordinate metadata can be recovered when both audited hg19 releases agree.
    coordinate_recovery<-semx_gwas_coordinates(ref,raws);ref<-coordinate_recovery$data
    data.table::fwrite(coordinate_recovery$audit,"reference_GWAS_coordinate_audit.csv")
    coordinate_mismatch<-ref$reason=="reference_coordinate_mismatch"
    ref$reason<-"included"
    invalid<-!is.finite(ref$MAF)|ref$MAF<.01|ref$MAF>.5|
      !ref$A1%in%c("A","C","G","T")|!ref$A2%in%c("A","C","G","T")|ref$A1==ref$A2
    ref$reason[invalid]<-"reference_MAF_or_allele_QC"
    bad_coordinates<-!ref$CHR%in%1:22|!is.finite(ref$BP)|ref$BP<=0
    ref$reason[!invalid&bad_coordinates]<-"reference_coordinates_unavailable_or_invalid"
    ref$reason[coordinate_mismatch]<-"reference_coordinate_mismatch"
    data.table::fwrite(ref,"reference_target_audit.csv")
    ref<-ref[ref$reason=="included",,drop=FALSE]
    writeLines(c(paste("Old reference rows:",old_reference_n),paste("Selected unique rsIDs:",length(ids)),
      paste("Selected rsIDs in old reference:",sum(ids%in%old_ids)),paste("Added reference rsIDs:",sum(!ref$old_reference)),
      paste("Coordinates recovered by concordant MDD/ANX:",sum(coordinate_recovery$audit$coordinate_action=="filled_from_concordant_MDD_ANX_GRCh37")),
      "The reference is audited empirically; its contents are not assumed to equal HapMap3."),"reference_summary.txt")
    for(trait in names(files)) {
      x<-raws[[trait]]
      pr<-semx_prepare_trait(x,trait,ref,ids);preps[[trait]]<-pr
      data.table::fwrite(pr$audit,paste0(trait,"_filter_audit.csv"))
      writeLines(capture.output(dput(pr$spec)),paste0(trait,"_column_mapping.txt"))
      data.table::fwrite(pr$data,paste0(trait,"_instruments_prepared.txt"),sep="\t",quote=FALSE,na="NA")
    }
    # Only independently frequency-resolved palindromes reach sumstats.
    common<-Reduce(intersect,lapply(preps,function(v) v$audit$SNP[v$audit$reason=="included"]))
    ref_ok<-ref[ref$SNP%in%common,,drop=FALSE]
    if(!nrow(ref_ok)) stop("No instrument SNP passed three-trait harmonization.")
    data.table::fwrite(ref_ok[,c("SNP","CHR","BP","MAF","A1","A2")],"reference_instruments_expanded.txt",sep="\t")
    q$ref_sign<-NA_real_;q$alignment<-"invalid_selected_pQTL"
    qef<-semx_col(names(q),c("EAF","EFFECT_ALLELE_FREQUENCY","A1FREQ"),FALSE)
    for(i in seq_len(nrow(q))) {
      if(!q$valid_pqtl[i]) next
      r<-match(q$SNP[i],ref_ok$SNP)
      if(is.na(r)) {q$alignment[i]<-"failed_three_trait_or_reference_QC";next}
      ef<-if(!is.na(qef)) semx_num(q[[qef]][i]) else NA_real_
      a<-semx_align(q$effect_allele[i],q$other_allele[i],ef,ref_ok$A1[r],ref_ok$A2[r],ref_ok$ref_EAF[r])
      q$alignment[i]<-a$action;if(a$ok) q$ref_sign[i]<-a$sign
    }
    q$eligible<-is.finite(q$ref_sign)
    data.table::fwrite(q,"pQTL_harmonization_audit.csv")
    all_audit<-Reduce(function(a,b) merge(a,b,by="SNP",all=TRUE),lapply(names(preps),function(t) {
      u<-preps[[t]]$audit[,c("SNP","raw_present","reason","raw_EAF","EAF_column"),drop=FALSE]
      names(u)[-1]<-paste0(t,"_",names(u)[-1]);u
    }))
    all_audit$in_old_reference<-all_audit$SNP%in%old_ids
    all_audit$in_expanded_QC_reference<-all_audit$SNP%in%ref_ok$SNP
    data.table::fwrite(all_audit,"all_instrument_filter_audit.csv")
    candidates<-c("AKT3","CD40","CGREF1","FES","FURIN","SIRPA","CD101","KHK","SNX18","UBE2L6")
    data.table::fwrite(merge(q[q$protein%in%candidates,,drop=FALSE],all_audit,by="SNP",all.x=TRUE),"candidate_filter_audit.csv")
    cat("Three-trait SNPs passing preparation=",nrow(ref_ok),"; frequency/allele-aligned proteins=",sum(q$eligible),"\n",sep="")
    if(isTRUE(getOption("ptsd.sem.expand.audit_only",FALSE))) {
      writeLines("Audit and input preparation complete; no new sumstats/userGWAS estimates.","AUDIT_COMPLETE.txt")
      return(invisible(out))
    }

    cat("STEP 4/6: Official sumstats on the expanded, vetted instrument reference\n")
    traits<-c("PTSD","MDD","ANX")
    call<-list(files=paste0(traits,"_instruments_prepared.txt"),ref="reference_instruments_expanded.txt",
      trait.names=traits,se.logit=rep(TRUE,3),OLS=rep(FALSE,3),linprob=rep(FALSE,3),
      N=c(NA_real_,829250,390020),info.filter=.6,maf.filter=.01,parallel=FALSE)
    # In the official implementation ambig=TRUE removes ambiguous SNPs.
    # Retention here is allowed only after the explicit EAF checks above.
    if("ambig"%in%names(formals(GenomicSEM::sumstats))) call$ambig<-FALSE
    ss<-as.data.frame(do.call(GenomicSEM::sumstats,call))
    if(anyDuplicated(ss$SNP)||!all(c("SNP","A1","A2")%in%names(ss))) stop("Unexpected sumstats output.")
    data.table::fwrite(ss,"expanded_instrument_sumstats.csv");saveRDS(ss,"expanded_instrument_sumstats.rds")
    ss_cols<-as.vector(rbind(paste0("beta.",traits),paste0("se.",traits)))
    if(!all(ss_cols%in%names(ss))) stop("Unexpected beta/SE names in sumstats output.")
    q$in_expanded_sumstats<-q$SNP%in%ss$SNP
    all_audit$in_expanded_sumstats<-all_audit$SNP%in%ss$SNP
    data.table::fwrite(all_audit,"all_instrument_filter_audit.csv")
    new_modes<-list()
    for(mode in c("primary","NEFF")) {
      cache<-base[[mode]]
      if(is.null(cache)) next
      od<-file.path(out,mode);dir.create(od)
      old_ss<-data.table::fread(file.path(cache$dir,"instrument_sumstats.csv"),data.table=FALSE)
      old_fg<-data.table::fread(file.path(cache$dir,"userGWAS_raw.csv"),data.table=FALSE)
      for(f in c("instrument_sumstats.csv","userGWAS_raw.csv","stage12_factorMR_results.csv"))
        file.copy(file.path(cache$dir,f),file.path(od,paste0("previous_",f)))
      previous<-cache$results
      previous_ids<-unique(previous$outcome_snp)
      if(!all(ss_cols%in%names(old_ss))) stop("Cached sumstats schema differs in ",mode)
      baseline_ok<-semx_same_ss(ss,old_ss,previous_ids,ss_cols)
      cmp<-data.frame(SNP=previous_ids,prepared_effects_and_SEs_match=baseline_ok)
      for(f in c("A1","A2","MAF",ss_cols)) {
        cmp[[paste0(f,"_previous")]]<-old_ss[[f]][match(previous_ids,old_ss$SNP)]
        cmp[[paste0(f,"_expanded")]]<-ss[[f]][match(previous_ids,ss$SNP)]
      }
      data.table::fwrite(cmp,file.path(od,"baseline_preparation_comparison.csv"))
      if(!all(baseline_ok)) stop(mode,": some of the original 391 proteins were lost or re-scaled; inspect baseline_preparation_comparison.csv before interpreting expanded results.")
      ld<-readRDS(file.path(cache$dir,"LDSCoutput_stage0.rds"))
      if(!is.matrix(ld$S)) ld$S<-as.matrix(ld$S)
      if(!all(dim(ld$S)==c(3L,3L))||!all(is.finite(ld$S))||!all(dim(ld$V)==c(6L,6L))) stop("Invalid cached S/V.")
      if(!is.null(colnames(ld$S))&&!identical(colnames(ld$S),traits)) stop("Cached trait order differs.")
      model<-paste(readLines(file.path(cache$dir,"model_primary.txt"),warn=FALSE),collapse="\n")
      expected<-"F1 =~ NA*PTSD + MDD + ANX\nF1 ~~ 1*F1\nPTSD ~~ PTSD\nMDD ~~ 0*MDD\nANX ~~ ANX"
      if(gsub("[[:space:]]","",model)!=gsub("[[:space:]]","",expected)) stop("Cached factor model differs from the reviewed constrained model.")
      fit_snps<-function(u,tag) {
        cat("  ",mode," ",tag,": ",nrow(u)," SNPs\n",sep="");flush.console()
        g<-GenomicSEM::userGWAS(covstruc=ld,SNPs=as.data.frame(u),estimation="DWLS",
          model=paste(model,"F1 ~ SNP",sep="\n"),sub="F1~SNP",Q_SNP=TRUE,cores=1,parallel=FALSE,printwarn=TRUE)
        saveRDS(g,file.path(od,paste0(tag,"_userGWAS.rds")))
        v<-as.data.frame(g[[1L]])
        if(!all(c("SNP","est","SE","Q_SNP","Q_SNP_pval")%in%names(v))||anyDuplicated(v$SNP)) stop("Unexpected userGWAS output.")
        data.table::fwrite(v,file.path(od,paste0(tag,"_userGWAS.csv")));v
      }
      cat("STEP 5/6: ",mode," control fits, then additional instrument SNPs\n",sep="")
      recovered_ids<-coordinate_recovery$audit$SNP[
        coordinate_recovery$audit$coordinate_action=="filled_from_concordant_MDD_ANX_GRCh37"]
      controls<-unique(c(head(previous_ids,6L),intersect(previous_ids,recovered_ids)))
      ctrl<-fit_snps(ss[match(controls,ss$SNP),,drop=FALSE],"control")
      ci<-match(ctrl$SNP,old_fg$SNP)
      control_ok<-!anyNA(ci)&&all(vapply(c("est","SE","Q_SNP","Q_SNP_pval"),function(f)
        all(is.finite(ctrl[[f]]))&&all(abs(ctrl[[f]]-old_fg[[f]][ci])<=1e-6*pmax(1,abs(old_fg[[f]][ci]))),logical(1)))
      writeLines(paste("Control model agreement:",control_ok),file.path(od,"control_agreement.txt"))
      if(!control_ok) stop(mode,": current userGWAS does not reproduce the cached control results. No mixing of model versions performed.")
      can_reuse<-semx_same_ss(ss,old_ss,ss$SNP,ss_cols)&ss$SNP%in%old_fg$SNP
      reused<-old_fg[match(ss$SNP[can_reuse],old_fg$SNP),,drop=FALSE]
      # Legacy userGWAS labels its first six input fields as SNP/CHR/BP/MAF/A1/A2.
      # The old four-column reference therefore produced malformed metadata only.
      # Repair these annotation fields from the current canonical sumstats;
      # cached coefficient, SE and Q values are preserved and checked above.
      for(f in c("SNP","CHR","BP","MAF","A1","A2"))
        reused[[f]]<-ss[[f]][match(reused$SNP,ss$SNP)]
      extra<-ss[!can_reuse,,drop=FALSE]
      added<-list()
      if(nrow(extra)) for(st in seq.int(1L,nrow(extra),by=100L)) {
        ii<-st:min(st+99L,nrow(extra));tag<-sprintf("added_%04d_%04d",min(ii),max(ii))
        added[[length(added)+1L]]<-fit_snps(extra[ii,,drop=FALSE],tag)
      }
      fg<-as.data.frame(data.table::rbindlist(c(list(reused),added),use.names=TRUE,fill=TRUE))
      good<-is.finite(fg$est)&is.finite(fg$SE)&fg$SE>0&is.finite(fg$Q_SNP)&fg$Q_SNP>=-1e-6&
        is.finite(fg$Q_SNP_pval)&fg$Q_SNP_pval>=0&fg$Q_SNP_pval<=1
      if("error"%in%names(fg)) good<-good&(is.na(fg$error)|as.character(fg$error)%in%c("0","FALSE",""))
      data.table::fwrite(fg,file.path(od,"all_userGWAS_raw.csv"))
      data.table::fwrite(fg[!good,,drop=FALSE],file.path(od,"failed_models.csv"))
      fg<-fg[good,,drop=FALSE]
      z<-q;ix<-match(z$SNP,fg$SNP)
      z$beta_factor<-fg$est[ix]*z$ref_sign;z$se_factor<-fg$SE[ix]
      z$Q<-fg$Q_SNP[ix];z$Q_pval<-fg$Q_SNP_pval[ix]
      z$testable<-z$eligible&z$in_expanded_sumstats&is.finite(z$beta_factor)&is.finite(z$se_factor)&z$se_factor>0
      z$reason<-ifelse(!z$eligible,z$alignment,ifelse(!z$in_expanded_sumstats,"excluded_by_sumstats",ifelse(!z$testable,"factor_model_unavailable_or_failed","included")))
      data.table::fwrite(z,file.path(od,"all_protein_audit.csv"))
      z<-z[z$testable,,drop=FALSE]
      z$mr_beta<-z$beta_factor/z$beta;z$mr_se<-abs(z$se_factor/z$beta);z$mr_z<-z$mr_beta/z$mr_se
      z$mr_p<-2*pnorm(-abs(z$mr_z));z$mr_fdr<-p.adjust(z$mr_p,"BH");z$Q_fdr<-p.adjust(z$Q_pval,"BH")
      z$outcome_snp<-z$SNP;z$exposure_EA<-z$effect_allele;z$exposure_OA<-z$other_allele
      z$pqtl_beta<-z$beta;z$pqtl_se<-z$se
      u<-!duplicated(z$SNP);ufdr<-p.adjust(z$Q_pval[u],"BH")
      z$Q_fdr_unique_SNP<-ufdr[match(z$SNP,z$SNP[u])]
      pi<-match(previous$protein,z$protein)
      if(anyNA(pi)||any(abs(z$mr_beta[pi]-previous$mr_beta)>1e-7)||any(abs(z$Q_pval[pi]-previous$Q_pval)>1e-6))
        stop(mode,": original protein estimates did not remain reproducible; inspect outputs.")
      comparison<-merge(previous[,c("protein","outcome_snp","mr_beta","mr_p","mr_fdr","Q_pval","Q_fdr")],
        z[,c("protein","outcome_snp","mr_beta","mr_p","mr_fdr","Q_pval","Q_fdr")],
        by=c("protein","outcome_snp"),all=TRUE,suffixes=c("_previous","_expanded"))
      data.table::fwrite(comparison,file.path(od,"comparison_with_previous.csv"))
      z<-z[order(z$mr_p),,drop=FALSE]
      data.table::fwrite(z,file.path(od,"factor_MR_expanded.csv"))
      data.table::fwrite(z[z$protein%in%candidates,,drop=FALSE],file.path(od,"candidate_results.csv"))
      counts<-data.frame(mode=mode,previous_proteins=nrow(previous),testable_proteins=nrow(z),unique_SNPs=length(unique(z$SNP)),
        MR_FDR_005=sum(z$mr_fdr<.05),Q_nominal_005=sum(z$Q_pval<.05),Q_FDR_005=sum(z$Q_fdr<.05),
        testable_candidates=sum(z$protein%in%candidates),newly_fitted_SNPs=nrow(extra),failed_models=sum(!good))
      data.table::fwrite(counts,file.path(od,"summary_counts.csv"));print(counts)
      new_modes[[mode]]<-counts
      dput(list(baseline_directory=cache$dir,ldsc_md5=unname(tools::md5sum(file.path(cache$dir,"LDSCoutput_stage0.rds"))),
        original_parameters=cache$pars,scope="all originally selected pQTL instruments passing documented QC",MAF=.01,INFO=.6,
        palindrome_max_MAF=.42,frequency_tolerance=.10,orientation_margin=.10,
        EAF_column_convention="Explicit A1 EAF, or PGC FRQ_U_<N> control frequency before FRQ_A_<N>; original headers preserved for package scaling.",
        coordinate_fallback="Concordant MDD/ANX GRCh37 positions and allele pairs when EUR coordinates unavailable.",
        note="Recorded LDSC/sample-size specifications reused. SNP variance, alleles and all trait effects/SEs matched the cached inputs for all original proteins. Malformed legacy userGWAS metadata repaired only; coefficient/SE/Q preserved. BH recomputed across the full expanded protein family."),
        file=file.path(od,"parameters.txt"))
    }
    cat("STEP 6/6: Save combined summaries\n")
    data.table::fwrite(data.table::rbindlist(new_modes,use.names=TRUE),"summary_counts_all_modes.csv")
    if(is.null(base$NEFF)) writeLines("NEFF LDSC cache unavailable: expanded source-NEFF sensitivity has not been run.","NEFF_NOT_RUN.txt")
    data.table::fwrite(merge(q[q$protein%in%candidates,,drop=FALSE],all_audit,by="SNP",all.x=TRUE),"candidate_filter_audit.csv")
    writeLines(c("COMPUTATION COMPLETE; inspect QC and control comparisons before changing manuscript counts.",
      "Pairwise Cochran Q is not Genomic SEM Q_SNP. Nominal heterogeneity alone is not corrected evidence of disorder specificity.",
      "Expand the multiplicity family consistently; preserve previous and expanded analyses as separate recorded versions."),"SUCCESS.txt")
    cat("Send candidate_filter_audit.csv, summary_counts_all_modes.csv, RUN.log, and\n",
        "primary/candidate_results.csv + primary/comparison_with_previous.csv.\n",sep="")
  },error=function(e) {
    writeLines(c(format(Sys.time()),conditionMessage(e)),file.path(out,"FAILED.txt"))
    cat("STOPPED: ",conditionMessage(e),"\n",sep="");stop(e)
  })
  invisible(out)
}
if(!isTRUE(getOption("ptsd.sem.expand.load_only",FALSE))) semx_run()
