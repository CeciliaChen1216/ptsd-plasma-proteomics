# Reconstruct the manuscript's Table S24c from the small, frozen audit outputs.
# Base R only. This does not read full GWAS data or run MR/LDSC/SEM/COJO.
# Cohort matching uses case/control deficits first; NEFF is checked afterwards.
# All count matches, including ambiguous ones, are retained.

ptsd_neff_pattern_data <- function(header_file, pattern_file, max_cohorts=3L) {
  if (length(max_cohorts)!=1L || !max_cohorts %in% 1:3)
    stop("max_cohorts must be 1, 2 or 3")
  lines <- readLines(header_file, warn=FALSE)
  field <- function(key) {
    prefix <- paste0("##", key, "=\"")
    hit <- lines[startsWith(lines, prefix)]
    if (length(hit)!=1L || !endsWith(hit, '"')) stop("Missing/ambiguous header field: ", key)
    substring(hit, nchar(prefix)+1L, nchar(hit)-1L)
  }
  cohorts <- strsplit(field("cohortList"), ",", fixed=TRUE)[[1]]
  cases <- as.numeric(strsplit(field("nCaseByCohort"), ",", fixed=TRUE)[[1]])
  controls <- as.numeric(strsplit(field("nControlByCohort"), ",", fixed=TRUE)[[1]])
  n_cohort <- as.integer(field("nCohort"))
  if (anyDuplicated(cohorts) || length(cohorts)!=n_cohort ||
      length(cases)!=n_cohort || length(controls)!=n_cohort ||
      any(!is.finite(cases)|cases<=0|!is.finite(controls)|controls<=0))
    stop("Invalid or inconsistent cohort header")
  header_cases <- as.numeric(field("nCase"))
  header_controls <- as.numeric(field("nControl"))
  if (sum(cases)!=header_cases || sum(controls)!=header_controls)
    stop("Cohort case/control sums differ from the release header")
  cohort <- data.frame(cohort=cohorts, cases=cases, controls=controls,
    binary_NEFF=4*cases*controls/(cases+controls), stringsAsFactors=FALSE)
  cohort$rounded_up_NEFF <- ceiling(cohort$binary_NEFF)

  patterns <- read.csv(pattern_file, check.names=FALSE, stringsAsFactors=FALSE)
  cols <- c("NCASE","NCON","NEFF","NTOT","SNP_rows")
  if (!all(cols %in% names(patterns)) || !nrow(patterns)) stop("Invalid count-pattern file")
  if (any(!is.finite(as.matrix(patterns[cols]))) ||
      any(as.matrix(patterns[cols])<=0) || anyDuplicated(patterns[cols[1:4]]))
    stop("Missing, non-positive or duplicated count patterns")
  if (nrow(patterns)>100L) stop("Expected at most the top 100 source patterns")
  # Preserve the original rank order, including ties.
  if (any(diff(patterns$SNP_rows)>0)) stop("Patterns are not ordered by frequency")
  ref <- patterns[1L,,drop=FALSE]

  # Enumerate each distinct combination once, using case/control counts only.
  combos <- unlist(lapply(seq_len(max_cohorts), function(k)
    combn(seq_len(n_cohort), k, simplify=FALSE)), recursive=FALSE)
  combo_cases <- vapply(combos, function(ix) sum(cases[ix]), numeric(1))
  combo_controls <- vapply(combos, function(ix) sum(controls[ix]), numeric(1))
  key <- function(a,b) paste(format(a, scientific=FALSE, trim=TRUE),
                           format(b, scientific=FALSE, trim=TRUE), sep="|")
  lookup <- split(seq_along(combos), key(combo_cases, combo_controls))
  rows <- list()
  for (i in seq_len(nrow(patterns))) {
    dc <- ref$NCASE-patterns$NCASE[i]
    du <- ref$NCON-patterns$NCON[i]
    dn <- ref$NEFF-patterns$NEFF[i]
    if (dc<0 || du<0 || (dc==0 && du==0)) next
    indices <- lookup[[key(dc,du)]]
    if (!length(indices)) next
    for (j in indices) {
      ix <- combos[[j]]
      predicted <- sum(cohort$rounded_up_NEFF[ix])
      rows[[length(rows)+1L]] <- data.frame(
        pattern_rank=i, SNP_rows=patterns$SNP_rows[i],
        matched_cohorts=paste(cohorts[ix], collapse=" + "),
        case_deficit=dc, control_deficit=du,
        observed_NEFF_deficit=dn,
        calculated_rounded_NEFF_deficit=predicted,
        observed_minus_calculated=dn-predicted,
        count_match_combinations=length(indices), stringsAsFactors=FALSE)
    }
  }
  matches <- if (length(rows)) do.call(rbind,rows) else data.frame(
    pattern_rank=integer(), SNP_rows=numeric(), matched_cohorts=character(),
    case_deficit=numeric(), control_deficit=numeric(), observed_NEFF_deficit=numeric(),
    calculated_rounded_NEFF_deficit=numeric(), observed_minus_calculated=numeric(),
    count_match_combinations=integer(), stringsAsFactors=FALSE)
  ranks <- unique(matches$pattern_rank)
  used_cohorts <- if (nrow(matches)) unique(unlist(strsplit(matches$matched_cohorts," + ",fixed=TRUE))) else character()
  exact_ranks <- unique(matches$pattern_rank[matches$observed_minus_calculated==0])
  list(cohorts=cohort, matches=matches, reference=ref,
    summary=list(header_cases=header_cases, header_controls=header_controls,
      n_header_cohorts=n_cohort, patterns_examined=nrow(patterns),
      top_pattern_SNP_rows=sum(patterns$SNP_rows),
      count_matched_patterns=length(ranks),
      count_matched_SNP_rows=sum(patterns$SNP_rows[ranks]),
      named_cohorts_in_matches=length(used_cohorts),
      exact_NEFF_matched_patterns=length(exact_ranks),
      ambiguous_count_patterns=length(unique(matches$pattern_rank[matches$count_match_combinations>1L])),
      header_binary_NEFF_sum=sum(cohort$binary_NEFF),
      header_rounded_NEFF_sum=sum(cohort$rounded_up_NEFF)))
}

ptsd_build_neff_tables <- function(audit_dir, output_dir, max_cohorts=3L) {
  needed <- c("source_header.txt","top100_source_N_patterns.csv",
    "source_count_consistency.csv","source_sample_size_summary.csv",
    "SE_approximation_agreement.csv","NEFF_formula_diagnostics.csv",
    "common_null_SNP_chromosome_diagnostics.csv","candidate_source_diagnostics.csv",
    "source_manifest.csv")
  missing <- needed[!file.exists(file.path(audit_dir,needed))]
  if (length(missing)) stop("Missing audit inputs: ",paste(missing,collapse=", "))
  manifest <- read.csv(file.path(audit_dir,"source_manifest.csv"),stringsAsFactors=FALSE)
  if (nrow(manifest)!=1L || !all(c("observed_md5","expected_md5") %in% names(manifest)) ||
      manifest$observed_md5!="bf0fb5b16bf91e640af14784e9d3a829" ||
      manifest$observed_md5!=manifest$expected_md5) stop("Unexpected source checksum record")
  read <- function(f) read.csv(file.path(audit_dir,f),check.names=FALSE,stringsAsFactors=FALSE)
  count <- read("source_count_consistency.csv")
  sample <- read("source_sample_size_summary.csv")
  agree <- read("SE_approximation_agreement.csv")
  formula <- read("NEFF_formula_diagnostics.csv")
  chrom <- read("common_null_SNP_chromosome_diagnostics.csv")
  candidate <- read("candidate_source_diagnostics.csv")
  one <- function(x,what) {if(length(x)!=1L || is.na(x)) stop("Missing/ambiguous audit value: ",what);x}
  check <- function(k) one(count$count[count$check==k],k)
  subset <- "common_MAF_ge_005_absZ_le_2"
  common <- agree[agree$subset==subset,,drop=FALSE]
  ratio <- formula[formula$subset==subset & formula$quantity=="SE_over_source_NEFF",,drop=FALSE]
  if (nrow(common)!=1L || nrow(ratio)!=1L) stop("Missing common-null diagnostics")
  value <- function(q,stat) one(sample[sample$quantity==q,stat],paste(q,stat))
  pat <- ptsd_neff_pattern_data(file.path(audit_dir,"source_header.txt"),
    file.path(audit_dir,"top100_source_N_patterns.csv"),max_cohorts)
  s <- pat$summary; ref <- pat$reference
  missing_checks <- c(check("NCASE_NCON_NTOT_missing_or_nonpositive"),check("NEFF_missing_or_nonpositive"))
  # A union count is identified by the separate totals only when both are zero.
  missing_union <- if(all(missing_checks==0)) 0 else "Not reconstructed from separate totals"
  labels <- c("Variants in the verified source",
    "Missing/non-positive NCASE, NCON, NTOT or NEFF",
    "NTOT vs NCASE + NCON discrepancies >1","Variants with NEFF > NTOT",
    "Header cases","Header controls","Source NTOT: minimum","Source NTOT: median",
    "Source NTOT: maximum","Source NEFF: minimum","Source NEFF: median","Source NEFF: maximum",
    "Common-null variants","Median N_SE / source NEFF",
    "N_SE / NEFF within 10% (common-null)","N_SE / NEFF within 20% (common-null)",
    "Smallest chromosome median N_SE / NEFF","Largest chromosome median N_SE / NEFF",
    "Prioritized instruments within 10%","Variants covered by the top 100 count patterns",
    "Patterns matching 1–3 header cohorts","Patterns with exact rounded-up NEFF deficit",
    "Variants in the 67 matched patterns","Header cohorts implicated by count matching",
    "Reference pattern NCASE","Reference pattern NCON","Reference pattern NEFF",
    "Unique global participant count established","All absolute cohort weights reconstructed")
  # Labels refer to the audited publication source; report an alternate search explicitly.
  if(max_cohorts!=3L) labels[21] <- paste0("Patterns matching 1–",max_cohorts," header cohorts")
  if(s$count_matched_patterns!=67L) labels[23] <- "Variants in matched count patterns"
  # In the audited source, max(NEFF) < min(NTOT); the audit also recorded a +1-tolerance check.
  neff_above_total <- if(value("NEFF","max")<=value("NTOT","min")) 0 else
    "Exact comparison not reconstructed; see original audit tolerance"
  vals <- list(check("source_rows"),missing_union,
    check("NTOT_differs_from_NCASE_plus_NCON_by_more_than_1"),neff_above_total,
    s$header_cases,s$header_controls,value("NTOT","min"),value("NTOT","median"),value("NTOT","max"),
    value("NEFF","min"),value("NEFF","median"),value("NEFF","max"),common$n,ratio$median,
    common$fraction_SE_ratio_within_010,common$fraction_SE_ratio_within_020,
    min(chrom$median_SE_over_source_NEFF),max(chrom$median_SE_over_source_NEFF),
    sum(is.finite(candidate$SE_over_source_NEFF) & abs(candidate$SE_over_source_NEFF-1)<=.1),
    s$top_pattern_SNP_rows,s$count_matched_patterns,s$exact_NEFF_matched_patterns,
    s$count_matched_SNP_rows,s$named_cohorts_in_matches,ref$NCASE,ref$NCON,ref$NEFF,"No","No")
  notes <- c("All source records","Count of affected variants",
    "Absolute difference in participants","Count of affected variants",
    "Release-level header; unique total unresolved","Release-level header; unique total unresolved",
    "Per-variant count, not unique global N","Per-variant count","Per-variant count, not unique global N",
    rep("Provider-reported information size",3),"MAF >=0.05 and |BETA/SE| <=2",
    "N_SE = 4/[2f(1-f)SE^2]",rep("Fraction of common-null variants",2),
    rep("Among 22 autosomes, common-null subset",2),"All ten prioritized instruments",
    "Count-pattern diagnostic, not unique participants","Exact NCASE and NCON deficits",
    "NEFF tested after case/control matching","Distinct source rows represented by these patterns",
    "Pattern matches do not prove cohort participation",
    "Most frequent pattern (1,986,989 variants)","Most frequent pattern",
    "Deficits measured against this pattern","Header/per-variant discrepancy remains unresolved",
    "Deficit matching is incomplete reconstruction")
  to_text <- function(x) if(is.numeric(x)) format(x,digits=17,scientific=FALSE,trim=TRUE) else x
  summary <- data.frame(Metric=labels,Value=vapply(vals,to_text,character(1)),
    Scope_and_interpretation=notes,stringsAsFactors=FALSE)
  if (!dir.exists(output_dir) && !dir.create(output_dir,recursive=TRUE)) stop("Cannot create output directory")
  files <- c("TableS24c_summary.csv","TableS24c_header_cohorts.csv",
    "TableS24c_count_pattern_matches.csv","TableS24c_pattern_diagnostics.csv","TableS24c_README.txt")
  if(any(file.exists(file.path(output_dir,files)))) stop("Refusing to overwrite existing S24c outputs")
  write <- function(x,f) write.csv(x,file.path(output_dir,f),row.names=FALSE,na="NA")
  write(summary,files[1]);write(pat$cohorts,files[2]);write(pat$matches,files[3])
  write(data.frame(metric=names(s),value=unlist(s,use.names=FALSE)),files[4])
  writeLines(c("Table S24c, Additional file 3: regenerated from frozen full-file audit outputs.",
    "Summary = workbook section A; header cohorts = B; count-pattern matches = C.",
    "The original R audit creates the small input files from the verified Figshare v3 GWAS.",
    "Search: one to three distinct header cohorts; exact case/control deficits relative to the most frequent pattern.",
    "All count matches are retained. NEFF is checked afterwards, not used to select matches.",
    "Cohort effective N = 4*cases*controls/(cases+controls), rounded upward separately within cohort.",
    "Count-pattern matching supports a cohort-additive structure but does not establish per-variant cohort membership.",
    "Absolute cohort weights and unique global participant counts remain unresolved.",
    "This export does not rerun MR/LDSC/SEM/COJO or rewrite the formatted Word/XLSX files."),
    file.path(output_dir,files[5]))
  invisible(list(summary=summary,cohorts=pat$cohorts,matches=pat$matches,diagnostics=s,
                 files=file.path(output_dir,files)))
}
