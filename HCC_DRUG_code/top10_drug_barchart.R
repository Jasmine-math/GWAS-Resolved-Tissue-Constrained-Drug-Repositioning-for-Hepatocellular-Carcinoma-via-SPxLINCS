# ---------- auto install packages ----------
pkgs <- c("readr", "dplyr", "stringr", "ggplot2", "forcats", "tidytext", "scales", "tibble")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(forcats)
library(tidytext)
library(scales)
library(grid)

# ---------- path list: ORIGINAL Table2 files ----------
paths <- c(
  "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/HCC_GCST90041897_hg38/Table2_selected_repositioning_hits_49tissues_HCC_parallel_HCC_GCST90041897_hg38.csv",
  "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/HCC_GCST90043858_hg38/Table2_selected_repositioning_hits_49tissues_HCC_parallel_HCC_GCST90043858_hg38.csv",
  "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/HCC_bbj-a-158_hg38/Table2_selected_repositioning_hits_49tissues_HCC_parallel_HCC_bbj-a-158_hg38.csv",
  "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/ieu-b-4915_hg38/Table2_selected_repositioning_hits_49tissues_HCC_parallel_ieu-b-4915_hg38.csv",
  "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/ebi-a-GCST90018583_hg38/Table2_selected_repositioning_hits_49tissues_HCC_parallel_ebi-a-GCST90018583_hg38.csv",
  "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/ebi-a-GCST90018803_hg38/Table2_selected_repositioning_hits_49tissues_HCC_parallel_ebi-a-GCST90018803_hg38.csv",
  "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/finn-b-CD2_BENIGN_LIVE_BILE_EXALLC_hg38/Table2_selected_repositioning_hits_49tissues_HCC_parallel_finn-b-CD2_BENIGN_LIVE_BILE_EXALLC_hg38.csv",
  "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/finn-b-CD2_BENIGN_LIVER_EXALLC_hg38/Table2_selected_repositioning_hits_49tissues_HCC_parallel_finn-b-CD2_BENIGN_LIVER_EXALLC_hg38.csv",
  "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/finn-b-C3_LIVER_INTRAHEPATIC_BILE_DUCTS_EXALLC_hg38/Table2_selected_repositioning_hits_49tissues_HCC_parallel_finn-b-C3_LIVER_INTRAHEPATIC_BILE_DUCTS_EXALLC_hg38.csv",
  "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/ebi-a-GCST90092003_hg38/Table2_selected_repositioning_hits_49tissues_HCC_parallel_ebi-a-GCST90092003_hg38.csv",
  "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/ebi-a-GCST90018858_hg38/Table2_selected_repositioning_hits_49tissues_HCC_parallel_ebi-a-GCST90018858_hg38.csv",
  "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/ebi-a-GCST90018638_hg38/Table2_selected_repositioning_hits_49tissues_HCC_parallel_ebi-a-GCST90018638_hg38.csv"
)

# ---------- trait labels ----------
trait_map <- c(
  "GCST90041897" = "Self reported hepatocellular carcinoma (Europe, UK Biobank)",
  "GCST90043858" = "ICD10 liver cell carcinoma (Europe, UK Biobank)",
  "bbj-a-158"    = "Hepatocellular carcinoma (Japan, BioBank)",
  "ieu-b-4915"   = "Liver and bile duct cancer (Europe)",
  "ebi-a-GCST90018583" = "Hepatic bile duct cancer (East Asia)",
  "ebi-a-GCST90018803" = "Hepatic bile duct cancer (Europe)",
  "finn-b-CD2_BENIGN_LIVE_BILE_EXALLC" = "Benign liver and bile ducts exallc (Europe, FinnGen)",
  "finn-b-CD2_BENIGN_LIVER_EXALLC"     = "Benign liver exallc (Europe, FinnGen)",
  "finn-b-C3_LIVER_INTRAHEPATIC_BILE_DUCTS_EXALLC" = "Liver intrahepatic bile ducts exallc (Europe, FinnGen)",
  "ebi-a-GCST90092003" = "Alcohol related hepatocellular carcinoma (Europe)",
  "ebi-a-GCST90018858" = "Hepatic cancer (Europe)",
  "ebi-a-GCST90018638" = "Hepatic cancer (East Asia)"
)

mapping_info <- tibble::tribble(
  ~id, ~group_row, ~hex,
  "GCST90041897", "R1", "#4c72b0",
  "GCST90043858", "R1", "#4c72b0",
  "bbj-a-158", "R1", "#4c72b0",
  "ieu-b-4915", "R2", "#dd8452",
  "ebi-a-GCST90018583", "R2", "#dd8452",
  "ebi-a-GCST90018803", "R2", "#dd8452",
  "finn-b-CD2_BENIGN_LIVE_BILE_EXALLC", "R3", "#55a868",
  "finn-b-CD2_BENIGN_LIVER_EXALLC", "R3", "#55a868",
  "finn-b-C3_LIVER_INTRAHEPATIC_BILE_DUCTS_EXALLC", "R3", "#55a868",
  "ebi-a-GCST90092003", "R4", "#c44e52",
  "ebi-a-GCST90018858", "R4", "#c44e52",
  "ebi-a-GCST90018638", "R4", "#c44e52"
)

trait_order <- stringr::str_wrap(as.character(trait_map[c(
  "GCST90041897", "GCST90043858", "bbj-a-158",
  "ieu-b-4915", "ebi-a-GCST90018583", "ebi-a-GCST90018803",
  "finn-b-CD2_BENIGN_LIVE_BILE_EXALLC", "finn-b-CD2_BENIGN_LIVER_EXALLC", "finn-b-C3_LIVER_INTRAHEPATIC_BILE_DUCTS_EXALLC",
  "ebi-a-GCST90092003", "ebi-a-GCST90018858", "ebi-a-GCST90018638"
)]), width = 35)

# ---------- utility functions ----------
canonical <- function(x) tolower(gsub("[\\s_\\-]+", "", x))

parse_id <- function(path) {
  fname <- basename(path)
  res <- sub("^.*_parallel_(.+?)_hg38\\.csv$", "\\1", fname)
  res <- sub("^HCC_", "", res)
  res
}

read_top10 <- function(path) {
  id <- parse_id(path)
  trait <- if (id %in% names(trait_map)) trait_map[[id]] else id
  
  message("Reading: ", path)
  
  df <- if (grepl("\t", readLines(path, n = 1))) {
    read_tsv(path, show_col_types = FALSE)
  } else {
    read_csv(path, show_col_types = FALSE)
  }
  
  names(df) <- gsub("^\ufeff", "", names(df))
  names(df) <- trimws(names(df))
  canon <- canonical(names(df))
  
  # find Drug column
  drug_idx <- which(canon == "drug")[1]
  if (is.na(drug_idx)) drug_idx <- 1
  
  # find q column
  q_idx <- which(canon %in% c("q", "qvalue", "qval", "fdr"))[1]
  if (is.na(q_idx)) q_idx <- which(grepl("q", canon) & grepl("value", canon))[1]
  if (is.na(q_idx)) stop("No q column found in: ", path)
  
  drug_col <- names(df)[drug_idx]
  q_col <- names(df)[q_idx]
  
  out <- df %>%
    mutate(
      q_num = suppressWarnings(as.numeric(.data[[q_col]])),
      Drug = as.character(.data[[drug_col]])
    ) %>%
    filter(!is.na(Drug), Drug != "") %>%
    group_by(Drug) %>%
    summarise(q_min = min(q_num, na.rm = TRUE), .groups = "drop") %>%
    filter(is.finite(q_min)) %>%
    arrange(q_min, Drug) %>%
    slice_head(n = 10) %>%
    mutate(
      q_min = ifelse(q_min <= 0, .Machine$double.xmin, q_min),
      neglog10_qmin = -log10(q_min),
      traitname = trait,
      gwas_id = id
    )
  
  return(out)
}

# ---------- read all ----------
df_all <- lapply(paths, read_top10) %>% bind_rows()

# ---------- prepare for plotting ----------
df_all <- df_all %>%
  mutate(gwas_id = trimws(gwas_id)) %>%
  left_join(mapping_info, by = c("gwas_id" = "id")) %>%
  mutate(
    traitname = stringr::str_wrap(traitname, width = 35),
    traitname = factor(traitname, levels = trait_order)
  ) %>%
  group_by(traitname) %>%
  mutate(drug_in_facet = tidytext::reorder_within(Drug, neglog10_qmin, traitname)) %>%
  ungroup()

if (any(is.na(df_all$hex))) {
  warning("Unmatched IDs for color mapping: ",
          paste(unique(df_all$gwas_id[is.na(df_all$hex)]), collapse = ", "))
}

# ---------- plot ----------
p_grid <- ggplot(df_all, aes(x = drug_in_facet, y = neglog10_qmin, fill = hex)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.1f", neglog10_qmin)),
            hjust = -0.2, size = 4.5, fontface = "bold") +
  coord_flip() +
  tidytext::scale_x_reordered() +
  facet_wrap(~ traitname, scales = "free_y", ncol = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  scale_fill_identity() +
  labs(
    title = "Drug Repurposing Candidates for HCC and Related Traits",
    x = "Candidate Compounds",
    y = expression(-log[10](q[min]))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    panel.spacing.x = unit(1.2, "lines"),
    panel.spacing.y = unit(1.2, "lines"),
    panel.grid.major.x = element_line(color = "grey90"),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18, margin = margin(b = 15))
  )

# ---------- save ----------
ggsave(
  "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/HCC_Drug_Grid_from_original_Table2_top10.png",
  p_grid, width = 18, height = 18, dpi = 300, bg = "white"
)

ggsave(
  "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/HCC_Drug_Grid_from_original_Table2_top10.pdf",
  p_grid, width = 18, height = 18, dpi = 300, bg = "white"
)