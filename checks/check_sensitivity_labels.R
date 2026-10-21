# Checks that reviewed classifications follow the displayed source statistics.
(function() {
 P <- getOption('ptsd.current.paths');if(is.null(P)) stop('Source 00_START.R first')
 rd <- function(n) read.csv(file.path(P$snapshot_dir,'sensitivity_review_20261019',n),stringsAsFactors=FALSE,na.strings=c('NA',''))
 e <- rd('eqtl_smr_direction.csv');r <- rd('reverse_mr_review.csv');m <- rd('evidence_matrix.csv')
 stopifnot(nrow(e)==20L,nrow(r)==10L,nrow(m)==80L)
 qualified <- is.finite(e$smr_fdr)&e$smr_fdr<.05&is.finite(e$p_HEIDI)&e$p_HEIDI>.05&is.finite(e$n_heidi)&e$n_heidi>=3
 stopifnot(identical(as.logical(e$qualified),qualified),sum(qualified)==9L,
   sum(qualified & e$direction=='Opposite')==5L,
   setequal(e$gene[qualified & e$direction=='Opposite'],c('KHK','SNX18','UBE2L6','FURIN')),
   sum(is.finite(e$p_SMR)&e$tissue=='brain')==8L,
   sum(is.finite(e$p_SMR)&e$tissue=='blood')==6L)
 tested <- is.finite(e$b_SMR)
 expect_dir <- ifelse(e$b_SMR[tested]*e$plasma_MR_beta[tested]>0,'Concordant','Opposite')
 stopifnot(identical(e$direction[tested],expect_dir))
 cls <- ifelse(r$IVW_P<.05 & r$WM_P<.05,'Both',ifelse(r$IVW_P<.05,'IVW only','Neither'))
 stopifnot(identical(r$reverse_class,cls),sum(r$IVW_P<.05)==6L,
   sum(cls=='Both')==1L,r$protein[cls=='Both']=='CD40',sum(cls=='IVW only')==5L)
 z <- m[m$domain=='Reverse-MR assessment',];z <- z[match(r$protein,z$protein),]
 stopifnot(identical(z$state,ifelse(cls=='Neither','I','C')))
 z <- m[m$domain=='Multi-instrument IVW',]
 stopifnot(z$state[z$protein=='SIRPA']=='I',z$state[z$protein=='CD101']=='C',sum(z$state=='N')==0L)
 z <- m[m$domain=='eQTL-SMR direction',]
 stopifnot(setequal(z$protein[z$state=='C'],c('KHK','SNX18','UBE2L6','FURIN')),
   setequal(z$protein[z$state=='S'],c('CD40','SIRPA','FES')))
 cat('PASS: transcript qualification/signs, reverse-estimator categories and tested multi-IVW labels agree with the reviewed source statistics.\n')
})()
