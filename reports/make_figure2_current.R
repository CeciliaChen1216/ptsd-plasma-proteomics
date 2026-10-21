# Plot-only wrapper: always use this release's reviewed expanded primary snapshot.
(function(){
 P <- getOption("ptsd.current.paths");if(is.null(P)) stop("Source 00_START.R first")
 old <- options(ptsd.figure2.input=file.path(P$snapshot_dir,"factor_MR.csv"),ptsd.figure2.output_dir=P$report_dir)
 on.exit(options(old),add=TRUE)
 source(file.path(P$package_root,"reports/plot_figure2_expanded.R"),local=new.env(parent=globalenv()),encoding="UTF-8")
})()
