P <- getOption("ptsd.current.paths")
if(is.null(P)) stop("Source 00_START.R first")
# PTSD source and provenance audit -- 2026-10-16
# Read-only audit of source data. No GWAS/MR/coloc/SEM analysis is run.
# Uses base R only. PDF text is included only if pdftools is already installed.
# The report records raw metadata and rows; it does not infer an effect scale.

ptsd_source_check <- function(
  input_dir = file.path(P$project_root,"data_source_check_20261016"),
  project_dir = P$project_root
) {
  if (!dir.exists(input_dir)) stop("Input directory does not exist: ", input_dir)
  out_base <- file.path(input_dir, paste0("audit_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  out <- out_base
  i <- 1L
  while (dir.exists(out)) {
    i <- i + 1L
    out <- paste0(out_base, "_", i)
  }
  if (!dir.create(out, recursive = TRUE)) stop("Cannot create output directory: ", out)
  report <- file.path(out, "PTSD_source_check_report.txt")
  report_con <- file(report, open = "wt", encoding = "UTF-8")
  sink_level <- sink.number()
  sink(report_con, split = TRUE)
  on.exit({
    while (sink.number() > sink_level) sink()
    close(report_con)
  }, add = TRUE)

  say <- function(...) cat(..., "\n", sep = "")
  section <- function(x) say("\n========== ", x, " ==========")
  safe_lines <- function(f) tryCatch(readLines(f, warn = FALSE),
    error = function(e) paste("READ ERROR:", conditionMessage(e)))
  show_file <- function(f) {
    if (!file.exists(f)) { say("NOT FOUND: ", f); return(invisible(NULL)) }
    say("FILE: ", f)
    say("BYTES: ", file.info(f)$size)
    say("MODIFIED: ", as.character(file.info(f)$mtime))
  }

  targets <- c(
    AKT3 = "rs4515770", CD101 = "rs12130298", CD40 = "rs4810485",
    CGREF1 = "rs2384572", FES = "rs1894401", FURIN = "rs2071410",
    KHK = "rs2304681", SIRPA = "rs6136377", SNX18 = "rs2548611",
    UBE2L6 = "rs28362950", prior_preview_1 = "rs641325",
    prior_preview_2 = "rs13387644"
  )
  # Delimiters prevent rs1894401 from matching rs18944010.
  make_pattern <- function(ids) paste0(
    "(^|[[:space:];,])(", paste(ids, collapse = "|"), ")([[:space:];,]|$)"
  )
  target_pattern <- make_pattern(unname(targets))

  scan_file <- function(f, label, expected_md5 = NULL) {
    section(label)
    show_file(f)
    if (!file.exists(f)) return(invisible(NULL))
    if (!is.null(expected_md5)) {
      say("Computing MD5 checksum...")
      observed_md5 <- unname(tools::md5sum(f))
      say("EXPECTED_MD5: ", expected_md5)
      say("OBSERVED_MD5: ", observed_md5)
      say("MATCHES_OFFICIAL_V3_BYTES: ", identical(tolower(observed_md5), tolower(expected_md5)))
      say("A checksum mismatch must be reviewed; it does not identify another version by itself.")
    }
    compressed <- grepl("\\.gz$", f, ignore.case = TRUE)
    if (compressed) {
      magic <- readBin(f, what = "raw", n = 2L)
      if (!identical(as.integer(magic), c(31L, 139L))) {
        say("ERROR: Not a gzip file. It may be an incomplete download or an HTML error page.")
        return(invisible(NULL))
      }
    }
    con <- if (compressed) gzfile(f, "rt") else file(f, "rt")
    on.exit(close(con), add = TRUE)
    headers <- character()
    first_rows <- character()
    matches <- character()
    issues <- character()
    n_body <- 0
    n_lines <- 0
    next_progress <- 2000000
    completed <- FALSE
    tryCatch({
      repeat {
        block <- withCallingHandlers(
          readLines(con, n = 100000L, warn = FALSE),
          warning = function(w) {
            issues <<- unique(c(issues, conditionMessage(w)))
            invokeRestart("muffleWarning")
          }
        )
        if (!length(block)) break
        n_lines <- n_lines + length(block)
        is_comment <- startsWith(block, "#")
        headers <- c(headers, block[is_comment])
        body <- block[!is_comment & nzchar(block)]
        if (length(first_rows) < 3L) first_rows <- head(c(first_rows, body), 3L)
        n_body <- n_body + length(body)
        keep <- grepl(target_pattern, body, perl = TRUE)
        if (any(keep)) matches <- c(matches, body[keep])
        if (n_body >= next_progress) {
          say("Scanned ", format(n_body, big.mark = ",", scientific = FALSE), " non-comment lines...")
          next_progress <- next_progress + 2000000
        }
      }
      completed <- TRUE
    }, error = function(e) {
      issues <<- unique(c(issues, paste("ERROR:", conditionMessage(e))))
    })
    say("FULL_STREAM_READ_COMPLETED: ", completed)
    say("TOTAL_LINES: ", format(n_lines, scientific = FALSE))
    say("NONEMPTY_NONCOMMENT_LINES: ", format(n_body, scientific = FALSE))
    say("For a flat file without a # header, the above count includes its column-name row.")
    say("\n--- METADATA / COLUMN HEADERS (verbatim) ---")
    if (length(headers)) say(paste(headers, collapse = "\n")) else say("No # comment/header lines.")
    say("\n--- FIRST THREE NON-COMMENT LINES (verbatim) ---")
    say(paste(first_rows, collapse = "\n"))
    say("\n--- KEY-SNP RECORDS (verbatim; duplicates retained) ---")
    if (length(matches)) say(paste(matches, collapse = "\n")) else say("No matching rsID tokens found.")
    say("\n--- MATCH COUNTS ---")
    for (j in seq_along(targets)) {
      n <- sum(grepl(make_pattern(targets[j]), matches, perl = TRUE))
      say(names(targets)[j], "\t", targets[j], "\t", n)
    }
    say("A missing rsID token is not proof that the genomic variant is absent.")
    say("Raw token matches require review of the ID column; no allele harmonization was applied.")
    if (length(issues)) {
      say("\nREAD WARNINGS / ERRORS:")
      say(paste(issues, collapse = "\n"))
    }
    invisible(NULL)
  }

  section("AUDIT DETAILS")
  say("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
  say("R: ", R.version.string)
  say("Input directory: ", input_dir)
  say("Output directory: ", out)
  say("Official release checked: https://figshare.com/articles/dataset/ptsd2024/26349322/3")
  say("Release DOI: 10.6084/m9.figshare.26349322.v3")
  say("Expected hashes below were obtained from the official Figshare API on 2026-10-16.")
  say("This report does not select a replacement GWAS, confirm log-odds units, or validate prior analyses.")

  section("OFFICIAL README")
  pdf <- file.path(input_dir, "ptsd2024_sumstats_readme.pdf")
  show_file(pdf)
  if (file.exists(pdf)) {
    pdf_md5 <- unname(tools::md5sum(pdf))
    say("OBSERVED_MD5: ", pdf_md5)
    say("EXPECTED_MD5: 7e666935eef19a5e7d1c71312311e791")
    say("MATCHES_OFFICIAL_V3_BYTES: ", identical(tolower(pdf_md5), "7e666935eef19a5e7d1c71312311e791"))
    if (requireNamespace("pdftools", quietly = TRUE)) {
      tryCatch({
        pages <- pdftools::pdf_text(pdf)
        for (j in seq_along(pages)) {
          say("\n--- README PAGE ", j, " ---")
          say(pages[j])
        }
      }, error = function(e) say("PDF_TEXT_ERROR: ", conditionMessage(e), "; upload the original PDF."))
    } else {
      say("PDF_TEXT_NOT_EXTRACTED: pdftools is not installed; no package was installed by this audit.")
      say("Please upload the original README PDF alongside this report.")
    }
  }

  scan_file(file.path(input_dir, "eur_ptsd_pcs_v4_aug3_2021.vcf.gz"),
    "NEW EUROPEAN OVERALL FILE", "033685ac31474f59e5a6ec297eca6b49")
  scan_file(file.path(input_dir, "eur_ptsdcasecontrol_pcs_v4_aug3_2021.vcf.gz"),
    "NEW EUROPEAN CASE-CONTROL FILE", "bf0fb5b16bf91e640af14784e9d3a829")

  section("PREVIOUSLY REPORTED RAW RECORDS -- USER CONSOLE OUTPUT")
  say("These are fingerprints from the previous console output, not newly verified source metadata.")
  say("ID\tA1\tA2\tFREQ\tBETA\tSE\tP\tNCASE\tNCON\tNEFF\tNTOT")
  say("rs1894401\tA\tG\t0.5351\t0.0243\t0.0043\t1.367e-08\t140667\t1108966\t460598\t1249633")
  say("rs2071410\tC\tG\t0.6894\t0.0177\t0.0049\t3.440e-04\t117906\t1067405\t404035\t1185311")
  say("Agreement at a few SNPs does not establish whole-file identity or ancestry.")
  scan_file(file.path(project_dir, "results", "ptsd_freeze3_smr.ma"),
    "OPTIONAL EXISTING CONVERTED .ma FILE")

  section("OPTIONAL EXISTING LD / COJO PROVENANCE")
  fams <- c(
    current_reference = file.path(P$workspace_root,"reference/1000G_v3/EUR.fam"),
    old_chr15_region = file.path(project_dir, "results", "cojo_chr15", "1000G_EUR_chr15_region.fam")
  )
  for (j in seq_along(fams)) {
    say("\n", names(fams)[j], ": ", fams[j])
    if (file.exists(fams[j])) {
      ff <- safe_lines(fams[j])
      if (any(startsWith(ff, "READ ERROR:"))) say(paste(ff, collapse = "\n")) else {
        say("NONEMPTY_FAM_ROWS: ", sum(nzchar(trimws(ff))))
        say("MD5: ", unname(tools::md5sum(fams[j])))
      }
    } else say("NOT FOUND")
  }
  log_dir <- file.path(project_dir, "results", "cojo_chr15")
  logs <- if (dir.exists(log_dir)) list.files(log_dir, pattern = "\\.log$", full.names = TRUE) else character()
  if (!length(logs)) say("No existing .log files found in: ", log_dir)
  for (lf in logs) {
    ll <- safe_lines(lf)
    kk <- which(grepl("bfile|keep|remove|individual|sample|subject|participant|cojo|GCTA|version|633|503|reference", ll, ignore.case = TRUE))
    say("\n--- ", basename(lf), " ---")
    if (length(kk)) {
      if (length(kk) > 100L) { say("Showing first 60 and last 40 matching log lines."); kk <- unique(c(head(kk, 60L), tail(kk, 40L))) }
      say(paste(paste0(kk, ": ", ll[kk]), collapse = "\n"))
    } else say("No selected provenance terms found.")
  }
  section("FINISHED")
  say("Audit finished. Inspect any NOT FOUND / WARNINGS / ERRORS above.")
  say("No source file was modified and no association analysis was rerun.")
  say("SEND THIS REPORT: ", report)
  say("ALSO SEND: ", pdf)
  invisible(list(output_dir = out, report = report))
}

if (!isTRUE(getOption("ptsd.source_check.load_only", FALSE))) ptsd_source_check()
