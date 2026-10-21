# Current package: path configuration only; scientific calculations preserved.
P <- getOption("ptsd.current.paths")
if (is.null(P)) stop("First source 00_START.R; then use ptsd_run().")
# PTSD Genomic SEM: isolated, uncached rerun. Created 2026-10-07.
# Run in a restarted R session: source(file.choose(), encoding = "UTF-8")
# Uses existing installed packages. Does not install/update packages.
# This is a reproducibility rerun under the existing study assumptions,
# not confirmation that prevalence/sample-size assumptions are optimal.

run_ptsd_genomicsem <- function() {
  root <- P$project_root
  out <- tempfile(pattern = paste0("genomicsem_rerun_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"),
                  tmpdir = file.path(root, "results"))
  if (!dir.exists(file.path(root, "results"))) stop("Missing D:/PTSD/results")
  if (!dir.create(out)) stop("Could not create new output folder")
  oldwd <- getwd()
  out_depth <- sink.number()
  msg_depth <- sink.number(type = "message")
  con <- file(file.path(out, "RUN.log"), open = "wt")
  sink(con, split = TRUE)
  sink(con, type = "message")
  on.exit({
    if (sink.number(type = "message") != msg_depth) sink(type = "message")
    while (sink.number() > out_depth) sink()
    close(con)
    setwd(oldwd)
    cat("\nOutput directory: ", out, "\n", sep = "")
  }, add = TRUE)
  tryCatch({
    cat("Start: ", format(Sys.time()), "\nOutput: ", out, "\n", sep = "")
    pkgs <- c("GenomicSEM", "data.table", "R.utils")
    missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing)) stop("Missing installed packages: ", paste(missing, collapse = ", "))
    suppressPackageStartupMessages(library(data.table))
    writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_start.txt"))
    writeLines(capture.output(packageDescription("GenomicSEM")), file.path(out, "GenomicSEM_description.txt"))
    writeLines(capture.output(dput(run_ptsd_genomicsem)), file.path(out, "rerun_function_snapshot.R.txt"))
    if (!"Q_SNP" %in% names(formals(GenomicSEM::userGWAS))) stop("Installed userGWAS does not expose Q_SNP. Send package description for review.")
    # Preserve the installed implementations, including development commit details.
    for (fn in c("munge", "ldsc", "sumstats", "usermodel", "userGWAS")) {
      writeLines(capture.output(dput(getExportedValue("GenomicSEM", fn))),
                 file.path(out, paste0("installed_", fn, ".R.txt")))
    }

    cfg <- new.env(parent = globalenv())
    config_file <- file.path(P$package_root, "00_config.R")
    if (file.exists(config_file)) sys.source(config_file, envir = cfg)
    pick <- function(key, fallback) {
      if (exists(key, envir = cfg, inherits = FALSE)) {
        x <- get(key, envir = cfg, inherits = FALSE)
        if (is.character(x) && length(x) == 1L && !grepl("/path/to/", x, fixed = TRUE)) return(x)
      }
      fallback
    }
    inputs <- c(
      PTSD = pick("ptsd_smr_ma", file.path(root, "results/ptsd_freeze3_smr.ma")),
      MDD = pick("mdd_path", file.path(root, "GWAS/PGC_MDD_2025/27061255/daner/daner_pgc_mdd_no23andMe-noUKBB_eur_hg19_v3.49.24.11.neff.gz")),
      ANX = pick("anx_path", file.path(root, "GWAS/PGC_Anxiety_2026/31389910/ANX_2026_daner_fullANX_v12_woUTAH_11022026.gz")),
      pQTL = pick("pqtl_path", file.path(root, "pQTL/ukbppp_cis_pqtl.rds")),
      HM3 = pick("hm3_path", P$hm3),
      reference = pick("ref_sumstats", P$sumstats_reference)
    )
    ld_dir <- pick("ldsc_dir", P$ldsc)
    ld_dir <- paste0(sub("[/\\\\]+$", "", ld_dir), "/")
    missing_files <- inputs[!file.exists(inputs)]
    if (length(missing_files)) stop("Input file(s) missing:\n", paste(names(missing_files), missing_files, collapse = "\n"))
    ld_files <- unlist(lapply(1:22, function(i) paste0(ld_dir, i, c(".l2.ldscore.gz", ".l2.M_5_50"))))
    if (any(!file.exists(ld_files))) stop("LD reference missing:\n", paste(ld_files[!file.exists(ld_files)], collapse = "\n"))
    # Enforce the intended phenotype releases; a different file needs review.
    if (basename(inputs["MDD"]) != "daner_pgc_mdd_no23andMe-noUKBB_eur_hg19_v3.49.24.11.neff.gz") stop("MDD input differs from the intended no23andMe-noUKBB release.")
    if (basename(inputs["ANX"]) != "ANX_2026_daner_fullANX_v12_woUTAH_11022026.gz") stop("ANX input differs from the intended woUTAH release.")
    setwd(out)
    manifest_paths <- c(unname(inputs), ld_files, if (file.exists(config_file)) config_file)
    cat("Hashing input files...\n")
    manifest <- data.frame(path = normalizePath(manifest_paths, winslash = "/"),
                           bytes = file.info(manifest_paths)$size,
                           md5 = unname(tools::md5sum(manifest_paths)))
    write.csv(manifest, "input_manifest.csv", row.names = FALSE)
    traits <- c("PTSD", "MDD", "ANX")
    sample_prev <- c(137136 / (137136 + 1085746), 0.5, 0.5)
    population_prev <- c(0.07, 0.15, 0.10)
    n_vec <- c(NA_real_, 829250, 390020)
    pars <- list(trait_order = traits, N = n_vec, sample.prev = sample_prev,
                 population.prev = population_prev, munge.info = 0.9,
                 sumstats.info = 0.6, MAF = 0.01, estimation = "DWLS",
                 note = "PTSD uses per-SNP total N; MDD and ANX use effective N. Existing study assumptions retained for this audit rerun.")
    dput(pars, file = "parameters.txt")
    print(pars)
    for (nm in c("MDD", "ANX")) {
      hdr <- fread(inputs[nm], nrows = 5)
      writeLines(capture.output(print(hdr)), paste0(nm, "_input_preview.txt"))
    }

    cat("\nSTEP 1/5: Prepare and munge source summary statistics\n")
    ptsd <- fread(inputs["PTSD"])
    need <- c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")
    if (!all(need %in% names(ptsd))) stop("PTSD .ma requires: ", paste(need, collapse = ", "))
    if (anyDuplicated(ptsd$SNP)) stop("Duplicate PTSD SNP IDs: resolve before running.")
    if (any(!is.finite(ptsd$N) | ptsd$N <= 0)) stop("Invalid PTSD per-SNP N")
    writeLines(capture.output(summary(ptsd$N)), "PTSD_N_summary.txt")
    setnames(ptsd, need, c("SNP", "A1", "A2", "MAF", "BETA", "SE", "P", "N"))
    prep <- file.path(out, "PTSD_prepared.txt")
    fwrite(ptsd, prep, sep = "\t")
    rm(ptsd); gc()
    gwas_files <- c(prep, unname(inputs[c("MDD", "ANX")]))
    GenomicSEM::munge(files = gwas_files, hm3 = unname(inputs["HM3"]),
                     trait.names = traits, N = n_vec, info.filter = 0.9, maf.filter = 0.01)
    munged <- paste0(traits, ".sumstats.gz")
    if (!all(file.exists(munged))) stop("Munge did not create all three files.")
    write.csv(data.frame(file = munged, md5 = unname(tools::md5sum(munged))),
              "munged_manifest.csv", row.names = FALSE)

    cat("\nSTEP 2/5: Fresh LDSC and trait models\n")
    ld <- GenomicSEM::ldsc(traits = munged, sample.prev = sample_prev,
                          population.prev = population_prev, ld = ld_dir, wld = ld_dir,
                          trait.names = traits)
    saveRDS(ld, "LDSCoutput_stage0.rds")
    S <- as.matrix(ld$S)
    write.csv(S, "LDSC_S.csv")
    write.csv(as.matrix(ld$V), "LDSC_V.csv")
    write.csv(cov2cor(S), "genetic_correlations.csv")
    if (!all(is.finite(S)) || any(diag(S) <= 0)) stop("Invalid LDSC covariance matrix.")
    # Compare the package's own liability-scale log against the returned matrix.
    logs <- list.files(out, pattern = "_ldsc\\.log$", full.names = TRUE)
    if (length(logs) == 1L) {
      lines <- readLines(logs, warn = FALSE)
      hl <- grep("Total Liability Scale h2:", lines, value = TRUE, fixed = TRUE)
      if (length(hl) == 3L) {
        h <- as.numeric(sub(".*Total Liability Scale h2: *([0-9.eE+-]+).*", "\\1", hl))
        check <- data.frame(trait = traits, log_h2 = h, matrix_h2 = diag(S), difference = diag(S) - h)
        write.csv(check, "LDSC_log_matrix_check.csv", row.names = FALSE)
        if (any(!is.finite(h)) || any(abs(h - diag(S)) > 0.00006)) stop("Fresh LDSC log and S matrix disagree; inspect LDSC_log_matrix_check.csv before continuing.")
      } else warning("Could not parse three liability h2 values; check the LDSC log manually.")
    } else warning("Could not identify a unique LDSC log; check manually.")
    model_free <- 'F1 =~ NA*PTSD + MDD + ANX
                   F1 ~~ 1*F1
                   PTSD ~~ PTSD
                   MDD ~~ MDD
                   ANX ~~ ANX'
    model_fixed <- 'F1 =~ NA*PTSD + MDD + ANX
                    F1 ~~ 1*F1
                    PTSD ~~ PTSD
                    MDD ~~ 0*MDD
                    ANX ~~ ANX'
    writeLines(model_fixed, "model_primary.txt")
    free_fit <- tryCatch(GenomicSEM::usermodel(covstruc = ld, estimation = "DWLS", model = model_free,
                         CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE),
                         error = function(e) list(error = conditionMessage(e)))
    fit <- GenomicSEM::usermodel(covstruc = ld, estimation = "DWLS", model = model_fixed,
                                    CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
    saveRDS(list(free = free_fit, constrained = fit), "stage0b_models.rds")
    writeLines(capture.output(print(fit)), "constrained_model_full.txt")
    writeLines(capture.output(print(free_fit)), "unconstrained_model_full.txt")
    if (is.data.frame(fit$results)) fwrite(fit$results, "constrained_model_parameters.csv")
    if (is.data.frame(fit$modelfit)) fwrite(fit$modelfit, "constrained_model_fit.csv")

    cat("\nSTEP 3/5: Fresh sumstats harmonization (may be slow)\n")
    ss <- GenomicSEM::sumstats(files = gwas_files, ref = unname(inputs["reference"]),
              trait.names = traits, se.logit = rep(TRUE, 3), OLS = rep(FALSE, 3),
              linprob = rep(FALSE, 3), N = n_vec, info.filter = 0.6, maf.filter = 0.01)
    saveRDS(ss, "stage1_sumstats.rds")
    ss <- as.data.table(ss)
    if (!all(c("SNP", "A1", "A2") %in% names(ss)) || anyDuplicated(ss$SNP)) stop("Invalid/duplicate harmonized SNP IDs")
    pqtl <- as.data.table(readRDS(inputs["pQTL"]))
    cols <- c("protein", "SNP", "beta", "se", "effect_allele", "other_allele")
    if (!all(cols %in% names(pqtl))) stop("Missing required pQTL columns: ", paste(setdiff(cols, names(pqtl)), collapse = ", "))
    if (anyDuplicated(pqtl$protein)) stop("More than one instrument per protein: resolve selection first.")
    ss_inst <- ss[SNP %in% pqtl$SNP]
    if (!nrow(ss_inst)) stop("No direct instruments in fresh sumstats")
    fwrite(ss_inst, "instrument_sumstats.csv")
    cat("Direct instrument SNPs: ", nrow(ss_inst), "\n", sep = "")

    cat("\nSTEP 4/5: Fresh instrument-only userGWAS; no cache\n")
    if (!"Q_SNP" %in% names(formals(GenomicSEM::userGWAS))) stop("Installed userGWAS does not expose Q_SNP. Send package description for review.")
    gw <- GenomicSEM::userGWAS(covstruc = ld, SNPs = as.data.frame(ss_inst),
          estimation = "DWLS", model = paste(model_fixed, "F1 ~ SNP", sep = "\n"),
          sub = "F1~SNP", Q_SNP = TRUE, cores = 1, parallel = FALSE, printwarn = TRUE)
    saveRDS(gw, "stage12_factor_snpeffects_userGWAS.rds")
    fg <- as.data.table(gw[[1]])
    fwrite(fg, "userGWAS_raw.csv")
    required <- c("SNP", "est", "SE", "Q_SNP", "Q_SNP_pval")
    if (!all(required %in% names(fg))) stop("Unexpected userGWAS columns; inspect userGWAS_raw.csv")
    if (anyDuplicated(fg$SNP)) stop("Duplicate userGWAS SNP rows")
    if ("error" %in% names(fg)) {
      bad <- !is.na(fg$error) & !as.character(fg$error) %in% c("0", "FALSE", "")
      if (any(bad)) stop("userGWAS contains model errors; inspect raw output")
    }
    if (any(!is.finite(fg$est) | !is.finite(fg$SE) | fg$SE <= 0 |
            !is.finite(fg$Q_SNP_pval) | fg$Q_SNP_pval < 0 | fg$Q_SNP_pval > 1)) stop("Invalid userGWAS estimates or Q P values")

    cat("\nSTEP 5/5: Full allele-pair matching, Wald MR and BH\n")
    key <- fg[, .(SNP, factor_beta_raw = est, se_factor = SE, Q = Q_SNP, Q_pval = Q_SNP_pval)]
    key <- merge(key, ss_inst[, .(SNP, ssA1 = toupper(A1), ssA2 = toupper(A2))], by = "SNP")
    a <- merge(pqtl, key, by = "SNP", all.x = TRUE)
    a[, `:=`(EA = toupper(effect_allele), OA = toupper(other_allele))]
    a[, direct := EA == ssA1 & OA == ssA2]
    a[, reverse := EA == ssA2 & OA == ssA1]
    a[, palindromic := paste0(EA, OA) %in% c("AT", "TA", "CG", "GC")]
    a[, reason := "included"]
    a[is.na(factor_beta_raw), reason := "No direct instrument in factor output"]
    a[reason == "included" & !(direct %in% TRUE | reverse %in% TRUE), reason := "Allele pair mismatch"]
    a[reason == "included" & (!is.finite(beta) | beta == 0 | !is.finite(se) | se <= 0), reason := "Invalid pQTL effect or SE"]
    # Palindromes are not resolved from allele labels alone. Exclude and report.
    a[reason == "included" & palindromic, reason := "Palindromic SNP: requires EAF/strand review"]
    fwrite(a, "instrument_audit.csv")
    fwrite(a[reason != "included", .(protein, SNP, reason)], "stage12_not_testable.csv")
    r <- a[reason == "included"]
    if (!nrow(r)) stop("No valid harmonized protein instruments")
    r[, beta_factor := factor_beta_raw * ifelse(direct, 1, -1)]
    r[, mr_beta := beta_factor / beta]
    r[, mr_se := abs(se_factor / beta)]
    r[, mr_z := mr_beta / mr_se]
    r[, mr_p := 2 * pnorm(-abs(mr_z))]
    r[, `:=`(mr_fdr = p.adjust(mr_p, "BH"), Q_fdr = p.adjust(Q_pval, "BH"))]
    if (!"cis_gene" %in% names(r)) r[, cis_gene := NA_character_]
    ans <- r[, .(protein, cis_gene, outcome_snp = SNP,
          exposure_EA = EA, exposure_OA = OA, sumstats_A1 = ssA1, sumstats_A2 = ssA2,
          pqtl_beta = beta, pqtl_se = se, beta_factor, se_factor,
          mr_beta, mr_se, mr_z, mr_p, mr_fdr, Q, Q_pval, Q_fdr)]
    setorder(ans, mr_p)
    fwrite(ans, "stage12_factorMR_results.csv")
    candidates <- c("AKT3", "CD40", "CGREF1", "FES", "FURIN", "SIRPA", "CD101", "KHK", "SNX18", "UBE2L6")
    fwrite(ans[protein %in% candidates], "candidate_results.csv")
    counts <- data.frame(testable_proteins = nrow(ans), unique_SNPs = uniqueN(ans$outcome_snp),
                        MR_FDR_005 = sum(ans$mr_fdr < .05), Q_nominal_005 = sum(ans$Q_pval < .05),
                        Q_FDR_005 = sum(ans$Q_fdr < .05))
    write.csv(counts, "summary_counts.csv", row.names = FALSE)
    print(counts)
    print(ans[protein %in% candidates])
    legacy_path <- file.path(root, "results/genomicsem/stage12_factorMR_results.csv")
    if (file.exists(legacy_path)) {
      legacy <- fread(legacy_path)
      cmp <- merge(ans[, .(protein, outcome_snp, mr_beta, mr_p, Q_pval)],
                   legacy[, .(protein, outcome_snp, mr_beta, mr_p, Q_pval)],
                   by = c("protein", "outcome_snp"), all = TRUE, suffixes = c("_rerun", "_legacy"))
      fwrite(cmp, "comparison_with_legacy.csv")
    }
    # Verify source inputs were not modified during the run.
    if (!identical(unname(tools::md5sum(manifest_paths)), manifest$md5)) stop("Input files changed during run; do not use results.")
    writeLines(capture.output(sessionInfo()), "sessionInfo_end.txt")
    writeLines(c("COMPUTATION COMPLETE; scientific review still required.", format(Sys.time())), "SUCCESS.txt")
    cat("\nCOMPUTATION COMPLETE. Upload summary_counts.csv, candidate_results.csv,\n",
        "LDSC_log_matrix_check.csv, constrained_model_full.txt and RUN.log for review.\n")
  }, error = function(e) {
    writeLines(c(format(Sys.time()), conditionMessage(e)), file.path(out, "FAILED.txt"))
    cat("\nSTOPPED: ", conditionMessage(e), "\n", sep = "")
    stop(e)
  })
  invisible(out)
}
run_ptsd_genomicsem()
