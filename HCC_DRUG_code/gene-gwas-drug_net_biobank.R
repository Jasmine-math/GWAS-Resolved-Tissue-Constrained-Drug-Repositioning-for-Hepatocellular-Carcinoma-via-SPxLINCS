# ====== Dependencies ======
options(repos = c(CRAN = "https://cran.rstudio.com"))
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(forcats)
})

# ====== Paths ======
sup6_path <- "D:/GWAS/HCC/spredixcan/eqtl/Supplementary Materials/Supplementary Tables 1-10/Supplementary Table 6.csv"

drug_files <- c(
  "HCC_bbj-a-158" =
    "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/HCC_bbj-a-158_hg38/Table2_selected_repositioning_hits_49tissues_HCC_parallel_HCC_bbj-a-158_hg38.csv",
  "HCC_GCST90041897" =
    "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/HCC_GCST90041897_hg38/Table2_selected_repositioning_hits_49tissues_HCC_parallel_HCC_GCST90041897_hg38.csv",
  "HCC_GCST90043858" =
    "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/HCC_GCST90043858_hg38/Table2_selected_repositioning_hits_49tissues_HCC_parallel_HCC_GCST90043858_hg38.csv"
)

trait_map <- c(
  "HCC_bbj-a-158"      = "Hepatocellular carcinoma (Japan, BioBank)",
  "HCC_GCST90041897"   = "Self reported hepatocellular carcinoma (Europe, UK Biobank)",
  "HCC_GCST90043858"   = "ICD10 liver cell carcinoma (Europe, UK Biobank)"
)

choose_traits <- unname(trait_map[names(drug_files)])  # Display names for the three GWAS

# ====== tools ======
read_auto <- function(path){
  first <- readLines(path, n=1)
  if (grepl("\t", first) && !grepl(",", first)) readr::read_tsv(path, show_col_types = FALSE)
  else readr::read_csv(path, show_col_types = FALSE)
}
canon <- function(x){
  x <- gsub("^\ufeff","", x); x <- trimws(x); tolower(gsub("[\\s_\\-\\.]+","", x))
}

# ====== 1) Gene Top-15 per trait from Supplementary Table 6 using nlogp ======
sup6 <- read_auto(sup6_path) %>%
  mutate(
    trait_label = trimws(trait_label),
    gene_name   = trimws(gene_name),
    nlogp       = suppressWarnings(as.numeric(nlogp)),
    pvalue      = suppressWarnings(as.numeric(pvalue)),
    zscore      = suppressWarnings(as.numeric(zscore))
  )

genes_top15_edges <- sup6 %>%
  filter(trait_label %in% choose_traits) %>%
  arrange(trait_label, desc(nlogp), pvalue) %>%
  group_by(trait_label, gene_name) %>%
  slice_max(order_by = nlogp, n = 1, with_ties = FALSE) %>%  # per trait, keep the record with max nlogp per gene
  ungroup() %>%
  group_by(trait_label) %>%
  slice_max(order_by = nlogp, n = 15, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    source_name = gene_name,     # Node label (Gene)
    target_trait = trait_label,  # Connects to GWAS (display name)
    edge_type   = "Gene->GWAS",
    nlogp, pvalue, zscore, tissue
  )

# ====== 2) Drug Top-15 per GWAS from full tables: min q + occurrence count ======
read_drug_top15 <- function(gwas_id, path, topN=15){
  df_raw <- read_auto(path)
  raw <- names(df_raw); can <- canon(raw)
  
  # Drug column
  drug_idx <- which(can %in% c("drug","drugname","compound","pertiname","name","pertiname"))
  if (length(drug_idx)==0) drug_idx <- 1
  drug_col <- raw[drug_idx[1]]
  
  # Prefer q; if missing, fall back to p
  q_idx <- which(can %in% c("q","qvalue","qval","fdr","padj","adjp","qmin"))
  p_idx <- which(can %in% c("p","pvalue","pval","p_val"))
  if (length(q_idx)>0) { val_col <- raw[q_idx[1]] } else if (length(p_idx)>0) { val_col <- raw[p_idx[1]] } else {
    # Fallback: any column containing q or p
    fb <- grep("q|p", can)
    if (length(fb)==0) stop("Cannot find q/p column: ", path)
    val_col <- raw[fb[1]]
  }
  
  df1 <- df_raw %>%
    mutate(
      Drug  =.data[[drug_col]],
      value = suppressWarnings(as.numeric(.data[[val_col]])),
      value = ifelse(!is.finite(value) | value <= 0, 1e-300, value)
    )
  
  #  Occurrence counts
  counts <- df1 %>% group_by(Drug) %>% summarise(count = n(),.groups = "drop")
  
  # Min value per drug; take Top-15
  top15 <- df1 %>%
    group_by(Drug) %>%
    summarise(q_min = min(value, na.rm = TRUE),.groups = "drop") %>%
    arrange(q_min) %>%
    distinct(Drug,.keep_all = TRUE) %>%
    slice_head(n = topN) %>%
    left_join(counts, by = "Drug") %>%
    mutate(
      neglog10_q  = -log10(q_min),
      trait_label = trait_map[[gwas_id]]
    ) %>%
    transmute(
      source_name = Drug,
      target_trait = trait_label,
      edge_type = "Drug->GWAS",
      q_min, neglog10_q, count
    )
  
  top15
}

drugs_top15_edges <- bind_rows(
  read_drug_top15("HCC_bbj-a-158",    drug_files[["HCC_bbj-a-158"]]),
  read_drug_top15("HCC_GCST90041897", drug_files[["HCC_GCST90041897"]]),
  read_drug_top15("HCC_GCST90043858", drug_files[["HCC_GCST90043858"]])
)

# ====== 3) Combine both edge types and build global nodes (shared gene/drug into unified nodes) ======
edges <- bind_rows(genes_top15_edges, drugs_top15_edges)

# Global node set: GWAS + unique genes + unique drugs
gwas_nodes <- tibble::tibble(
  shared_name = paste0("gwas:", choose_traits),
  label = choose_traits,
  type  = "GWAS"
)
gene_nodes <- edges %>%
  filter(edge_type=="Gene->GWAS") %>%
  distinct(source_name) %>%
  transmute(shared_name = paste0("gene:", source_name),
            label = source_name, type = "Gene")
drug_nodes <- edges %>%
  filter(edge_type=="Drug->GWAS") %>%
  distinct(source_name) %>%
  transmute(shared_name = paste0("drug:", source_name),
            label = source_name, type = "Drug")

nodes <- bind_rows(gwas_nodes, gene_nodes, drug_nodes)

# ====== 4) Compute three-layer coordinates (fixed variable names) ======
x_map <- setNames(c(-300, 0, 300), choose_traits)
y_gwas <- 0; y_gene <- 200; y_drug <- -200

# Set GWAS coordinates
nodes <- nodes %>%
  mutate(`X Location` = NA_real_, `Y Location` = NA_real_) %>%
  mutate(
    `X Location` = ifelse(type=="GWAS", unname(x_map[label]), `X Location`),
    `Y Location` = ifelse(type=="GWAS", y_gwas, `Y Location`)
  )

# Gene node coordinates (use target_trait)）
gene_x <- edges %>%
  filter(edge_type=="Gene->GWAS") %>%
  transmute(shared_name = paste0("gene:", source_name),
            target_trait) %>% 
  mutate(xg = unname(x_map[target_trait])) %>%
  group_by(shared_name) %>%
  summarise(x0 = mean(xg, na.rm = TRUE), .groups="drop") %>%
  arrange(x0) %>%
  group_by(grp = round(x0, 0)) %>%
  mutate(offset = if (n()==1) 0 else seq(-80, 80, length.out = n())) %>%
  ungroup() %>%
  transmute(shared_name, `X Location` = x0 + offset, `Y Location` = y_gene)

# Drug node coordinates (use target_trait)
drug_x <- edges %>%
  filter(edge_type=="Drug->GWAS") %>%
  transmute(shared_name = paste0("drug:", source_name),
            target_trait) %>%
  mutate(xg = unname(x_map[target_trait])) %>%
  group_by(shared_name) %>%
  summarise(x0 = mean(xg, na.rm = TRUE), .groups="drop") %>%
  arrange(x0) %>%
  group_by(grp = round(x0, 0)) %>%
  mutate(offset = if (n()==1) 0 else seq(-80, 80, length.out = n())) %>%
  ungroup() %>%
  transmute(shared_name, `X Location` = x0 + offset, `Y Location` = y_drug)

# Merge coordinates robustly to avoid rows_update issues
nodes_update <- bind_rows(gene_x, drug_x) %>%
  rename(new_x = `X Location`, new_y = `Y Location`)  # avoid suffix collisions

# Left join and coalesce to fill coordinates
nodes <- nodes %>%
  left_join(nodes_update, by = "shared_name") %>%
  mutate(
    `X Location` = coalesce(new_x, `X Location`),
    `Y Location` = coalesce(new_y, `Y Location`)
  ) %>%
  select(shared_name, label, type, `X Location`, `Y Location`)

# ====== 5) Build final edge table ======
edges_out <- edges %>%
  mutate(
    source = paste0(ifelse(edge_type=="Gene->GWAS", "gene:", "drug:"), source_name),
    target = paste0("gwas:", target_trait),
    width_metric = ifelse(edge_type=="Gene->GWAS", nlogp, neglog10_q),
    color_metric = ifelse(edge_type=="Gene->GWAS", nlogp, neglog10_q)
  ) %>%
  select(
    source, target, edge_type,
    nlogp, pvalue, zscore, tissue,
    q_min, neglog10_q, count,
    width_metric, color_metric
  )

# ====== 6) Export ======
out_nodes <- "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/cyto_tripartite_nodes_HCC3.csv"
out_edges <- "D:/Drug/HCC/hcc drug paper final data & code/SPxLINCS_results_strict7169/cyto_tripartite_edges_HCC3.csv"

write_csv(nodes, out_nodes)
write_csv(edges_out, out_edges)
cat("Saved:\n", out_nodes, "\n", out_edges, "\n")