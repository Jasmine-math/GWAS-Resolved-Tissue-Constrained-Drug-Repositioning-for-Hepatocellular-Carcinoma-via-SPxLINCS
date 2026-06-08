# ==============================================================================
# Publication-Ready SIDER Indication Enrichment Bubble Plot (English-only)
# ==============================================================================

library(ggplot2)
library(dplyr)

sider_path <- "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/SIDER_Indications_table.txt"

if (file.exists(sider_path)) {
  
  # Load and preprocess SIDER table
  raw_sider <- read.delim(sider_path, sep = "\t", stringsAsFactors = FALSE)
  raw_sider$Count <- as.numeric(gsub("/.*", "", raw_sider$Overlap))
  raw_sider$Combined.Score <- as.numeric(raw_sider$Combined.Score)
  
  # 1) Select top 15 most significant clinical indications by Combined Score
  plot_df <- raw_sider %>%
    arrange(desc(Combined.Score)) %>%
    head(15)
  
  # 2) Text cleanup for overly long “stage unspecified” terms for better readability
  plot_df$Term <- gsub("hodgkin's disease lymphocyte predominance type stage unspecified",
                       "Hodgkin's disease (lymphocyte predominance)", plot_df$Term, ignore.case = TRUE)
  plot_df$Term <- gsub("hodgkin's disease lymphocyte depletion type stage unspecified",
                       "Hodgkin's disease (lymphocyte depletion)", plot_df$Term, ignore.case = TRUE)
  
  # Capitalize first letter for a more formal medical appearance
  plot_df$Term <- paste0(toupper(substr(plot_df$Term, 1, 1)), substr(plot_df$Term, 2, nchar(plot_df$Term)))
  
  # Reorder factor by Combined Score for plotting
  plot_df$Term <- reorder(plot_df$Term, plot_df$Combined.Score)
  
  # 3) Build a publication-quality bubble plot
  p_sider_final <- ggplot(plot_df, aes(x = Combined.Score, y = Term)) +
    geom_segment(aes(x = 0, xend = Combined.Score, y = Term, yend = Term),
                 color = "gray90", size = 0.5) +
    geom_point(aes(size = Count, color = Combined.Score), alpha = 0.9) +
    scale_color_gradientn(colors = c("#41B6C4", "#225EA8", "#081D58")) +
    scale_size_continuous(range = c(5, 11)) +
    labs(
      x = "Combined Score (Drug-Set Enrichment Analysis)",
      y = "Enriched SIDER Clinical Indications",
      title = "Clinical Disease Indications Enriched by 103 Candidate Drugs",
      size = "Overlapped Drug Count",
      color = "Enrichment Magnitude"
    ) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.major.x = element_line(color = "gray95"),
      panel.grid.major.y = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      # Tuning text size and vertical justification to avoid label overlap
      axis.text.y = element_text(color = "black", face = "bold", size = 10, vjust = 0.5),
      legend.position = "right"
    )
  
  print(p_sider_final)
  
  # Save the final high-resolution figure
  output_png <- "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/HCC_True_DSEA_SIDER_Final.png"
  ggsave(output_png, plot = p_sider_final, width = 12, height = 8, dpi = 300, type = "cairo")
  cat("[Done] Final publication-quality SIDER bubble plot has been saved.\n")
  
} else {
  cat("[Error] Input file not found. Please check the 'sider_path'.\n")
}