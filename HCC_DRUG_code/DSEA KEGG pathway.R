# ==============================================================================
# 1. Environment
# ==============================================================================
library(ggplot2)
library(dplyr)

kegg_down_path <- "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/L1000FWD_KEGG_Pathways_Down_table.txt"
kegg_up_path   <- "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/L1000FWD_KEGG_Pathways_Up_table.txt"

# ==============================================================================
# 2. Load data with precise column alignment
# ==============================================================================

# 2.1 Process Down‑regulated pathways (to be shown on the left; negative score)
raw_down <- read.delim(kegg_down_path, sep = "\t", stringsAsFactors = FALSE) %>% head(10)

data_down <- data.frame(
  Pathway = raw_down$Term,                        # Pathway term
  Score   = -as.numeric(raw_down$Combined.Score), # Use negative Combined Score for "Down"
  Pvalue  = as.numeric(raw_down$P.value),         # P value
  Count   = raw_down$Overlap,                     # Overlap column (e.g., "5/120")
  Regulation = "Down",
  stringsAsFactors = FALSE
)

# 2.2 Process Up‑regulated pathways (to be shown on the right; positive score)
raw_up <- read.delim(kegg_up_path, sep = "\t", stringsAsFactors = FALSE) %>% head(10)

data_up <- data.frame(
  Pathway = raw_up$Term,
  Score   = as.numeric(raw_up$Combined.Score),    # Keep positive Combined Score for "Up"
  Pvalue  = as.numeric(raw_up$P.value),
  Count   = raw_up$Overlap,
  Regulation = "Up",
  stringsAsFactors = FALSE
)

# ==============================================================================
# 3. Tidy and clean
# ==============================================================================
plot_df <- rbind(data_down, data_up)

# Clean Count: extract the numerator from strings like "5/120"
plot_df$Count <- as.numeric(gsub("/.*", "", plot_df$Count))

# Factorize and lock the y‑axis order by Score
plot_df$Regulation <- factor(plot_df$Regulation, levels = c("Down", "Up"))
plot_df$Pathway <- reorder(plot_df$Pathway, plot_df$Score)

# ==============================================================================
# 4. Draw a bidirectional bubble scatter plot
# ==============================================================================
p <- ggplot(plot_df, aes(x = Score, y = Pathway)) +
  # Midline at zero
  geom_vline(xintercept = 0, color = "gray50", linetype = "dashed", size = 0.6) +
  # Bubbles: size = gene count, color = regulation direction
  geom_point(aes(size = Count, color = Regulation), alpha = 0.85) +
  # Academic color palette
  scale_color_manual(values = c("Up" = "#D62728", "Down" = "#1F77B4")) +
  # Bubble size range
  scale_size_continuous(range = c(5, 11)) +
  # Labels
  labs(
    x = "Combined Score\n(← Inhibited by Drugs | Activated by Drugs →)",
    y = "Enriched KEGG Pathways",
    title = "KEGG Pathways Regulated by 103 Candidate HCC Drugs",
    size = "Target Gene Count",
    color = "Drug Regulation"
  ) +
  # Publication-grade theming
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.x = element_line(color = "gray95"),
    panel.grid.major.y = element_line(color = "gray90", linetype = "dotted"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    axis.title = element_text(face = "bold"),
    axis.text.y = element_text(color = "black", face = "bold", size = 10),
    legend.position = "right",
    legend.box.background = element_rect(color = "gray80", size = 0.5)
  )

# Preview in R
print(p)

# ==============================================================================
# 5. Save high‑resolution figures
# ==============================================================================
output_pdf <- "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/HCC_Drugs_KEGG_Dual_BubblePlot.pdf"
ggsave(output_pdf, plot = p, width = 11, height = 7.5, dpi = 300)

output_png <- "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/HCC_Drugs_KEGG_Dual_BubblePlot.png"
ggsave(output_png, plot = p, width = 11, height = 7.5, dpi = 300, type = "cairo")

cat("[Done] High‑resolution KEGG dual bubble plot saved to:\n", output_pdf, "\n")