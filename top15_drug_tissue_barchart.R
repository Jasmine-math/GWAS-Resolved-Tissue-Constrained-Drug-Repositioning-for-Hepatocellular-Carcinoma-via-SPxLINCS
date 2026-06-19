suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(tidytext)   # reorder_within/scale_y_reordered
  library(scales)     # label_wrapped
})

# 1)12 GWAS catalogue（Table2 CSV）
dirs <- c(
  "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/finn-b-C3_LIVER_INTRAHEPATIC_BILE_DUCTS_EXALLC_hg38",
  "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/finn-b-CD2_BENIGN_LIVE_BILE_EXALLC_hg38",
  "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/ebi-a-GCST90018583_hg38",
  "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/finn-b-CD2_BENIGN_LIVER_EXALLC_hg38",
  "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/HCC_bbj-a-158_hg38",
  "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/HCC_GCST90041897_hg38",
  "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/HCC_GCST90043858_hg38",
  "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/ieu-b-4915_hg38",
  "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/ebi-a-GCST90018638_hg38",
  "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/ebi-a-GCST90018803_hg38",
  "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/ebi-a-GCST90018858_hg38",
  "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/ebi-a-GCST90092003_hg38"
)

# 2) find each Table2
pattern <- "^Table2_selected_repositioning_hits_49tissues_HCC_parallel_.*\\.csv$"
files <- unlist(lapply(dirs, function(d) list.files(d, pattern = pattern, full.names = TRUE)))
stopifnot(length(files) > 0)

# 3) read each Table2 and get GWAS_ID
read_one <- function(fp){
  gwas_id <- str_match(basename(fp),
                       "^Table2_selected_repositioning_hits_49tissues_HCC_parallel_(.+)\\.csv$")[,2]
  readr::read_csv(fp, show_col_types = FALSE) |>
    rename(
      Drug   = any_of(c("Drug")),
      Tissue = any_of(c("Tissue")),
      q      = any_of(c("q value","q_value","q")),
      P      = any_of(c("P value","P_value","P","p","p_value","p value"))
    ) |>
    mutate(GWAS_ID = ifelse(is.na(gwas_id), "NA_GWAS", gwas_id),
           nlog10q = -log10(q))
}

df_all <- purrr::map_dfr(files, read_one)

# 4) in each Tissue×GWAS×Drug keep min q and only remain one when repeated
df_dedup <- df_all |>
  group_by(Tissue, GWAS_ID, Drug) |>
  slice_min(q, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(nlog10q = -log10(q))

# 5) summarize through GWAS：in each Tissue×Drug get min q and calculate repeated times
agg_df <- df_dedup |>
  group_by(Tissue, Drug) |>
  summarise(
    q_min_across_gwas = min(q, na.rm = TRUE),
    nlog10q           = -log10(q_min_across_gwas),
    recurrence        = n_distinct(GWAS_ID),
    .groups = "drop"
  )

# 6) in each Tissue get Top-15 and save（49）
out_dir <- "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/combined/per_tissue_top15"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

make_plot_one_tissue <- function(tis, top_n = 15, width = 7){
  df_t <- agg_df |>
    filter(Tissue == tis) |>
    slice_min(q_min_across_gwas, n = top_n, with_ties = FALSE) |>
    arrange(q_min_across_gwas) |>
    mutate(Drug_label = gsub("_", " ", Drug)) 
  
  if (nrow(df_t) == 0L) return(invisible(NULL))

  h <- max(4, 0.4 * nrow(df_t) + 1.2)
  
  p <- ggplot(
    df_t,
    aes(x = nlog10q,
        y = tidytext::reorder_within(Drug_label, nlog10q, Tissue),
        fill = recurrence)
  ) +
    geom_col() +
    tidytext::scale_y_reordered(
      labels = function(x) {
        x <- sub("___.*$", "", x)         
        x <- stringr::str_squish(x)        
        scales::label_wrap(28)(x)          
      }
    ) +
    scale_fill_viridis_c(option = "C", name = "GWAS\nrecurrence") +
    labs(
      title = tis,
      x = expression(-log[10]("min q across GWAS")),
      y = "Drug (Top 15)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 10),
      legend.position = "right"
    )
  
  fn <- file.path(out_dir, paste0("Top15_", gsub("[^A-Za-z0-9_]+","_", tis), ".png"))
  ggsave(fn, p, width = width, height = h, dpi = 300)
  message("Saved: ", fn)
  fn <- file.path(out_dir, paste0("Top15_", gsub("[^A-Za-z0-9_]+","_", tis), ".pdf"))
  ggsave(fn, p, width = width, height = h, dpi = 300)
  message("Saved: ", fn)
}

tissues <- sort(unique(agg_df$Tissue))
invisible(lapply(tissues, make_plot_one_tissue))

message("All per-tissue figures written to: ", out_dir)