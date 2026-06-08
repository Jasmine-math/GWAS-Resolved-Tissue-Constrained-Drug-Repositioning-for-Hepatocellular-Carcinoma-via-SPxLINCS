# -*- coding: utf-8 -*-
import pandas as pd
from pathlib import Path
import numpy as np

# =========================
# 1. GWAS result directories
# =========================
base_dirs = [
    r"D:\Drug\HCC\hcc drug paper final data & code\SPxLINCS_results_strict7169\finn-b-CD2_BENIGN_LIVER_EXALLC_hg38",
    r"D:\Drug\HCC\hcc drug paper final data & code\SPxLINCS_results_strict7169\HCC_bbj-a-158_hg38",
    r"D:\Drug\HCC\hcc drug paper final data & code\SPxLINCS_results_strict7169\HCC_GCST90041897_hg38",
    r"D:\Drug\HCC\hcc drug paper final data & code\SPxLINCS_results_strict7169\HCC_GCST90043858_hg38",
    r"D:\Drug\HCC\hcc drug paper final data & code\SPxLINCS_results_strict7169\ieu-b-4915_hg38",
    r"D:\Drug\HCC\hcc drug paper final data & code\SPxLINCS_results_strict7169\ebi-a-GCST90018583_hg38",
    r"D:\Drug\HCC\hcc drug paper final data & code\SPxLINCS_results_strict7169\ebi-a-GCST90018638_hg38",
    r"D:\Drug\HCC\hcc drug paper final data & code\SPxLINCS_results_strict7169\ebi-a-GCST90018803_hg38",
    r"D:\Drug\HCC\hcc drug paper final data & code\SPxLINCS_results_strict7169\ebi-a-GCST90018858_hg38",
    r"D:\Drug\HCC\hcc drug paper final data & code\SPxLINCS_results_strict7169\ebi-a-GCST90092003_hg38",
    r"D:\Drug\HCC\hcc drug paper final data & code\SPxLINCS_results_strict7169\finn-b-C3_LIVER_INTRAHEPATIC_BILE_DUCTS_EXALLC_hg38",
    r"D:\Drug\HCC\hcc drug paper final data & code\SPxLINCS_results_strict7169\finn-b-CD2_BENIGN_LIVE_BILE_EXALLC_hg38",
]

# =========================
# 2. Key tissues
# =========================
KEY_TISSUES_STRONG = {"Liver"}
KEY_TISSUES_MID = {
    "Whole_Blood",
    "Cells_EBV-transformed_lymphocytes",
    "Spleen",
    "Lung",
    "Small_Intestine_Terminal_Ileum",
    "Colon_Sigmoid",
    "Colon_Transverse"
}
KEY_TISSUES = KEY_TISSUES_STRONG | KEY_TISSUES_MID

# =========================
# 3. Processing function
# =========================
def process_one_gwas(gwas_dir):
    gwas_dir = Path(gwas_dir)
    gwas_name = gwas_dir.name

    # Find Table2 file
    pattern = f"Table2_selected_repositioning_hits_49tissues_HCC_parallel_{gwas_name}.csv"
    inp_path = gwas_dir / pattern

    if not inp_path.exists():
        print(f"[SKIP] File not found: {inp_path}")
        return None, None

    print(f"\n[PROCESS] {gwas_name}")
    print(f"Input -> {inp_path}")

    # Read
    df = pd.read_csv(inp_path)

    # Required columns
    required_cols = ["Tissue", "q value", "P value", "Z_agg_tissue"]
    for c in required_cols:
        if c not in df.columns:
            print(f"[SKIP] Missing column '{c}' in {gwas_name}")
            return None, None

    # Numeric conversion
    for c in ["q value", "P value", "Z_agg_tissue"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # Filter key tissues
    df_key = df[df["Tissue"].isin(KEY_TISSUES)].copy()

    # Strict reverse hits
    df_key_rev = df_key[df_key["Z_agg_tissue"] < 0].copy()

    # Sort within tissue
    df_key_rev = df_key_rev.sort_values(
        ["Tissue", "q value", "P value", "Z_agg_tissue"],
        ascending=[True, True, True, True]
    ).copy()

    # Top15 per tissue
    top15_list = []
    q15_values = []

    for t, sub in df_key_rev.groupby("Tissue"):
        sub15 = sub.head(15).copy()
        sub15["selection_stage"] = "top15_per_tissue"
        top15_list.append(sub15)

        if len(sub15) >= 15 and pd.notnull(sub15.iloc[14]["q value"]):
            q15_values.append(sub15.iloc[14]["q value"])

    if top15_list:
        df_top15 = pd.concat(top15_list, axis=0, ignore_index=True)
    else:
        df_top15 = pd.DataFrame(columns=df_key_rev.columns)

    # Global threshold
    if len(q15_values) == 0:
        threshold_q = df_top15["q value"].min() if len(df_top15) > 0 else np.inf
        print(f"[WARN] {gwas_name}: no tissue has >=15 rows; threshold = {threshold_q}")
    else:
        threshold_q = min(q15_values)
        print(f"[INFO] {gwas_name}: threshold_q = {threshold_q}")

    # Apply threshold
    df_sel = df_top15[df_top15["q value"] <= threshold_q].copy()
    df_sel["selection_reason"] = "q<=min_q_of_15th"

    # Fallback top3 for tissues with no selected rows
    tissues_all = set(df_top15["Tissue"].unique())
    tissues_have = set(df_sel["Tissue"].unique())
    tissues_lack = sorted(list(tissues_all - tissues_have))

    fallback_rows = []
    for t in tissues_lack:
        sub_all = df_key_rev[df_key_rev["Tissue"] == t].copy()
        fallback = sub_all.head(3).copy()
        fallback["selection_reason"] = "fallback_top3"
        fallback["selection_stage"] = "fallback_from_all_reverse"
        fallback_rows.append(fallback)

    # Combine
    if fallback_rows:
        df_final = pd.concat([df_sel] + fallback_rows, axis=0, ignore_index=True)
    else:
        df_final = df_sel.copy()

    df_final = df_final.sort_values(["Tissue", "q value", "P value"], ascending=[True, True, True])

    # Keep useful columns
    cols = [
        "Drug", "Cell line", "Tissue", "Rank in tissue", "P value", "q value",
        "Brief description", "n_contexts", "Z_agg_tissue", "selection_stage", "selection_reason"
    ]
    cols = [c for c in cols if c in df_final.columns]
    df_final = df_final[cols]

    # Output file
    out_path = gwas_dir / f"KEY_TISSUES_STRONG_MID_strict_top15_then_threshold_fallback3_{gwas_name}.csv"
    df_final.to_csv(out_path, index=False)
    print(f"Saved -> {out_path}")

    # Summary
    if len(df_final) > 0 and "Drug" in df_final.columns:
        summary = (
            df_final.groupby("Tissue")["Drug"]
            .count()
            .rename("n_rows")
            .reset_index()
            .sort_values(["Tissue"], ascending=[True])
        )
    else:
        summary = pd.DataFrame(columns=["Tissue", "n_rows"])

    summary["gwas"] = gwas_name
    summary["threshold_q"] = threshold_q

    return df_final, summary

# =========================
# 4. Run all GWAS
# =========================
all_summaries = []

for d in base_dirs:
    df_final, summary = process_one_gwas(d)
    if summary is not None:
        all_summaries.append(summary)

# =========================
# 5. Save merged summary
# =========================
if all_summaries:
    merged_summary = pd.concat(all_summaries, axis=0, ignore_index=True)

    summary_out = Path(r"D:\Drug\HCC\hcc drug paper final data & code\SPxLINCS_results_strict7169") / \
                  "SUMMARY_key_tissues_topdrug_counts_across_12GWAS.csv"
    merged_summary.to_csv(summary_out, index=False)
    print(f"\n[ALL DONE] Summary saved -> {summary_out}")
else:
    print("\n[DONE] No summary generated.")