# Historical classification retained for provenance; use the reviewed current plot.
stop("This historical matrix is superseded. Source 00_START.R and run ptsd_reports(suppfig2=TRUE).")
# Path compatibility only; load the current entry point before this historical module.
P <- getOption("ptsd.current.paths")
if(is.null(P)) stop("Source 00_START.R first")
###############################################################################
# SuppFig2_evidence_matrix_STANDALONE.R
# Generates Supplementary Figure S2 only.
# No 00_config.R or external result CSV is required.
###############################################################################

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required. Install it with install.packages('ggplot2').")
}

library(ggplot2)

## Output folder ---------------------------------------------------------------
output_dir <- file.path(P$project_root,"Figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## Evidence states -------------------------------------------------------------
# S = Supportive; C = Cautionary; I = Inconclusive; N = Not testable
# Protein order matches the final Table S1.
protein_order <- c(
  "AKT3", "CGREF1", "FES", "KHK", "SNX18",
  "CD40", "FURIN", "UBE2L6", "SIRPA", "CD101"
)

axis_labels <- c(
  "Primary\nMR",
  "Steiger\ndirectionality",
  "No MR-PRESSO\npleiotropy",
  "No reverse-MR\nevidence",
  "Multi-instrument\nIVW",
  "eQTL-SMR+\nHEIDI",
  "Colocalization",
  "Brain pQTL/\nPWAS"
)

# Rows = proteins; columns = evidence domains.
M <- rbind(
  AKT3   = c("S", "S", "S", "S", "S", "N", "S", "C"),
  CGREF1 = c("S", "S", "S", "S", "S", "I", "S", "C"),
  FES    = c("S", "S", "S", "S", "S", "S", "N", "I"),
  KHK    = c("S", "S", "S", "S", "S", "S", "S", "S"),
  SNX18  = c("S", "S", "S", "S", "S", "S", "N", "I"),
  CD40   = c("S", "S", "S", "C", "S", "S", "N", "S"),
  FURIN  = c("S", "S", "S", "S", "S", "S", "N", "N"),
  UBE2L6 = c("S", "S", "S", "S", "S", "S", "C", "I"),
  SIRPA  = c("S", "S", "S", "S", "N", "S", "N", "C"),
  CD101  = c("S", "S", "C", "S", "N", "N", "N", "N")
)
colnames(M) <- axis_labels

state_names <- c(
  S = "Supportive",
  C = "Cautionary",
  I = "Inconclusive",
  N = "Not testable"
)

plot_data <- as.data.frame(as.table(M), stringsAsFactors = FALSE)
names(plot_data) <- c("Protein", "Evidence domain", "State")

plot_data$Protein <- factor(
  plot_data$Protein,
  levels = rev(protein_order)
)
plot_data$`Evidence domain` <- factor(
  plot_data$`Evidence domain`,
  levels = axis_labels
)
plot_data$State <- factor(
  unname(state_names[plot_data$State]),
  levels = c("Supportive", "Cautionary", "Inconclusive", "Not testable")
)

palette <- c(
  "Supportive"   = "#2C7FB8",
  "Cautionary"   = "#D95F0E",
  "Inconclusive" = "#BDBDBD",
  "Not testable" = "#F0F0F0"
)

## Plot ------------------------------------------------------------------------
p <- ggplot(
  plot_data,
  aes(x = `Evidence domain`, y = Protein, fill = State)
) +
  geom_tile(colour = "white", linewidth = 1.05) +
  scale_fill_manual(values = palette, drop = FALSE, name = NULL) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(
    title = "Evidence matrix for ten PTSD-prioritized proteins",
    subtitle = paste(
      "Eight complementary evidence domains;",
      "cells are interpreted separately and not summed"
    ),
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x.top = element_text(
      face = "bold",
      size = 8.4,
      lineheight = 0.95,
      hjust = 0.5,
      margin = margin(b = 5)
    ),
    axis.text.y = element_text(face = "bold", size = 10),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.key.width = grid::unit(1.25, "cm"),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 9, colour = "grey30", hjust = 0.5),
    plot.margin = margin(t = 12, r = 18, b = 8, l = 8)
  )

## Export ----------------------------------------------------------------------
ggsave(
  file.path(output_dir, "SuppFig2_evidence_matrix.pdf"),
  p, width = 11.5, height = 7.4
)
ggsave(
  file.path(output_dir, "SuppFig2_evidence_matrix.tiff"),
  p, width = 11.5, height = 7.4,
  dpi = 300, compression = "lzw"
)
ggsave(
  file.path(output_dir, "SuppFig2_evidence_matrix.png"),
  p, width = 11.5, height = 7.4,
  dpi = 300
)

write.csv(
  plot_data,
  file.path(output_dir, "SuppFig2_evidence_matrix_states.csv"),
  row.names = FALSE
)

cat(
  "Supplementary Figure S2 written to:\n",
  normalizePath(output_dir, winslash = "/", mustWork = FALSE),
  "\n"
)
