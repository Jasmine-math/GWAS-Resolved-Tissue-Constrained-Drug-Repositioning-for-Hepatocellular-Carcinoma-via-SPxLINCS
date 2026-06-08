# -*- coding: utf-8 -*-
"""
End-to-end SPxLINCS refinement (English version)
- Translated all comments, logs, and section headers to English
- Keeps original functionality and parameterization
"""
import os, re, time, argparse, warnings, sys
import numpy as np
import pandas as pd
from glob import glob
from pathlib import Path
from scipy.stats import spearmanr, pearsonr, norm
from statsmodels.stats.multitest import multipletests
import h5py
import multiprocessing as mp
from datetime import datetime

# Unbuffered stdout for real-time logs
os.environ.setdefault("PYTHONUNBUFFERED", "1")
try:
    sys.stdout.reconfigure(line_buffering=True)
except Exception:
    pass

warnings.filterwarnings("ignore", category=FutureWarning)

# ================= Key-tissue definitions =================
# Core set of 8 tissues (Liver + 7 immune/gastrointestinal tissues)
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
KEY_TISSUES_CORE = KEY_TISSUES_STRONG | KEY_TISSUES_MID  # total of 8

# ================= Utilities =================
def ts():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def winsorize(s, q=0.01):
    lo, hi = s.quantile(q), s.quantile(1 - q)
    return s.clip(lo, hi)

def zscore_series(s):
    s = s.astype(float)
    sd = s.std(ddof=0)
    return (s - s.mean()) / sd if sd > 0 else s * 0

def stouffer(pvals, signs=None, weights=None):
    pvals = np.clip(np.array(pvals, dtype=float), 1e-300, 1 - 1e-16)
    Z = norm.isf(pvals / 2.0)
    if signs is not None:
        Z = Z * np.sign(np.array(signs, dtype=float))
    if weights is None:
        weights = np.ones_like(Z)
    Z_agg = np.sum(weights * Z) / np.sqrt(np.sum(weights**2))
    p_agg = 2 * norm.sf(abs(Z_agg))
    return Z_agg, p_agg

def _read_table_auto(path):
    sep = "\t" if str(path).lower().endswith((".txt", ".tsv", ".txt.gz", ".tsv.gz")) else ","
    return pd.read_csv(path, sep=sep, compression="infer", low_memory=False)

def build_gctx_index(gctx_path):
    with h5py.File(gctx_path, "r") as f:
        row_ids = f["0"]["META"]["ROW"]["id"][...]
        col_ids = f["0"]["META"]["COL"]["id"][...]
        row_ids = [x.decode("utf-8") if hasattr(x, "decode") else str(x) for x in row_ids]
        col_ids = [x.decode("utf-8") if hasattr(x, "decode") else str(x) for x in col_ids]
        rid2idx = {rid: i for i, rid in enumerate(row_ids)}
        cid2idx = {cid: i for i, cid in enumerate(col_ids)}
    return row_ids, col_ids, rid2idx, cid2idx

def find_data_dataset(f):
    try:
        node = f["0"]["DATA"]["0"]
        if isinstance(node, h5py.Dataset):
            return node
        for k in getattr(node, "keys", lambda: [])():
            if isinstance(node[k], h5py.Dataset):
                return node[k]
    except Exception:
        pass
    holder = []
    def visitor(name, obj):
        if isinstance(obj, h5py.Dataset) and len(obj.shape) == 2:
            holder.append(obj)
    f.visititems(visitor)
    if not holder:
        raise RuntimeError("No 2D dataset found in GCTX.")
    return holder[0]

def load_gctx_cols_windowed(gctx_path, row_idx_list, col_idx_list, row_ids, rid2symbol_map, window_size=200, lock=None):
    row_idx = np.unique(np.array(row_idx_list, dtype=np.int64)); row_idx.sort()
    orig_pos = np.arange(len(col_idx_list))
    pairs = sorted(zip(orig_pos, col_idx_list), key=lambda x: x[1])
    sorted_pos = np.array([p[0] for p in pairs], dtype=np.int64)
    sorted_cols = np.array([p[1] for p in pairs], dtype=np.int64)

    windows = []
    for i in range(0, len(sorted_cols), window_size):
        win_cols = sorted_cols[i:i + window_size]
        win_pos = sorted_pos[i:i + window_size]
        windows.append((win_cols, win_pos))

    if lock is not None:
        lock.acquire()
    try:
        with h5py.File(gctx_path, "r") as f:
            dset = find_data_dataset(f)
            shape = dset.shape
            submat = np.empty((len(row_idx), len(col_idx_list)), dtype=np.float32)

            if len(row_ids) == shape[0]:
                for (win_cols, win_pos) in windows:
                    c0, c1 = int(win_cols.min()), int(win_cols.max())
                    block = dset[row_idx, slice(c0, c1 + 1)].astype(np.float32)
                    rel = (win_cols - c0).astype(np.int64)
                    submat[:, win_pos] = block[:, rel]
                row_rids = np.array(row_ids)[row_idx]
            elif len(row_ids) == shape[1]:
                for (win_cols, win_pos) in windows:
                    c0, c1 = int(win_cols.min()), int(win_cols.max())
                    block = dset[slice(c0, c1 + 1), row_idx].astype(np.float32).T
                    rel = (win_cols - c0).astype(np.int64)
                    submat[:, win_pos] = block[:, rel]
                row_rids = np.array(row_ids)[row_idx]
            else:
                raise RuntimeError(f"row_ids length ({len(row_ids)}) mismatches dataset shape {shape}")
    finally:
        if lock is not None:
            lock.release()

    row_symbols = [rid2symbol_map.get(rid, rid) for rid in row_rids]
    df = pd.DataFrame(submat, index=row_symbols)
    df = df.groupby(df.index).mean()
    return df

def select_gene_sets_from_Z(z_series, k):
    s = z_series.dropna().sort_values(ascending=False)
    k_eff = min(k, len(s) // 2)
    return set(s.index[:k_eff]), set(s.index[-k_eff:])

def _es_from_hit_indices(w, hit_idx):
    N = w.shape[0]
    hits = np.zeros(N, dtype=bool); hits[hit_idx] = True
    w_hit_sum = float(np.sum(w[hits])); n_miss = int((~hits).sum())
    RS = np.empty(N, dtype=float); ph = 0.0; pm = 0
    w_hit_den = (w_hit_sum if w_hit_sum > 0 else 1.0)
    miss_den = (n_miss if n_miss > 0 else 1.0)
    for i in range(N):
        if hits[i]:
            ph += w[i]
            RS[i] = (ph / w_hit_den) - (pm / miss_den)
        else:
            pm += 1
            RS[i] = (ph / w_hit_den) - (pm / miss_den)
    return RS.max() - RS.min()

def ks_connectivity_fast(z_series, e_series, k, alpha, nperm, rng=None):
    rng = rng or np.random.default_rng(42)
    genes = z_series.index.intersection(e_series.index)
    if len(genes) < max(50, k):
        return np.nan, np.nan, np.nan, len(genes)

    z = z_series.loc[genes].values
    e = e_series.loc[genes].values

    order_z = np.argsort(-z)
    k_eff = min(k, len(order_z) // 2)
    up_idx = order_z[:k_eff]
    dn_idx = order_z[-k_eff:]

    order_e = np.argsort(-e)
    w = np.abs(e[order_e]) ** alpha
    pos_in_e = np.empty(len(e), dtype=int); pos_in_e[order_e] = np.arange(len(e))

    up_hit = pos_in_e[up_idx]
    dn_hit = pos_in_e[dn_idx]

    ES_up = _es_from_hit_indices(w, np.sort(up_hit))
    ES_dn = _es_from_hit_indices(w, np.sort(dn_hit))
    CS = ES_up - ES_dn

    N = len(w); m_up = len(up_hit); m_dn = len(dn_hit)
    null = np.empty(nperm, dtype=float)
    for i in range(nperm):
        perm_up = rng.choice(N, size=m_up, replace=False)
        perm_dn = rng.choice(N, size=m_dn, replace=False)
        ES_up_p = _es_from_hit_indices(w, np.sort(perm_up))
        ES_dn_p = _es_from_hit_indices(w, np.sort(perm_dn))
        null[i] = ES_up_p - ES_dn_p

    mu0 = float(np.mean(null))
    sd0 = float(np.std(null, ddof=1)) if np.std(null, ddof=1) > 0 else 1.0
    tau = (CS - mu0) / sd0
    p_two = ((np.sum(np.abs(null - mu0) >= np.abs(CS - mu0)) + 1) / (nperm + 1)) * 2
    return tau, min(p_two, 1.0), CS, len(genes)

# ================= Stage-2 multiprocessing worker =================
def refine_one_tissue(args_tuple):
    (
        t, df_t, gctx, row_idx_list, row_ids, rid2symbol_map,
        k_fast, nperm_fast, top_pct, k_full, nperm_strong, nperm_mid,
        alpha, batch_size, winsor_q, min_genes, sig_info_dict, selected_cids, cid2idx,
        window_size, sp_by_tissue_dict, force_full, h5_lock
    ) = args_tuple

    t_start = time.time()

    # Mode selection:
    # - Core 8 tissues -> full (K=50,100,250,500)
    # - Other tissues -> fast (K=50,100); top_pct=1.0 means all go into Stage 2
    mode = "fast"; K_list = k_fast; nperm = nperm_fast
    if force_full:
        mode = "full"; K_list = k_full; nperm = nperm_strong
    else:
        if t in KEY_TISSUES_STRONG:
            mode = "full"; K_list = k_full; nperm = nperm_strong
        elif t in KEY_TISSUES_MID:
            mode = "full"; K_list = k_full; nperm = nperm_mid

    if df_t.empty:
        return t, [], f"[WARN] {t}: no screened pairs"

    if mode == "fast":
        df_t = df_t.sort_values("p_screen")
        n_sel = max(1, int(len(df_t) * top_pct))  # top_pct=1.0 -> all
        sel = df_t.iloc[:n_sel]
        non = df_t.iloc[n_sel:]
        sel_sig_ids = sel["sig_id"].unique().tolist()
        logmsg = f"[INFO] {t} FAST (ALL): sel={len(sel_sig_ids)}/{len(df_t)}, K={K_list}, nperm={nperm}"
    else:
        sel = df_t
        non = df_t.iloc[0:0]
        sel_sig_ids = sel["sig_id"].unique().tolist()
        logmsg = f"[INFO] {t} FULL: sel={len(sel_sig_ids)}, K={K_list}, nperm={nperm}"

    final_rows = []
    calc_cnt = 0

    for sub_chunk in (sel_sig_ids[i:i + batch_size] for i in range(0, len(sel_sig_ids), batch_size)):
        col_idx_list = np.array([cid2idx[cid] for cid in sub_chunk], dtype=np.int64)

        expr_df = load_gctx_cols_windowed(
            gctx, row_idx_list, col_idx_list, row_ids, rid2symbol_map, window_size=window_size, lock=h5_lock
        )
        expr_df = expr_df.apply(lambda s: zscore_series(winsorize(s, q=winsor_q)))

        for j, sig_id in enumerate(sub_chunk):
            mm = sig_info_dict.get(sig_id, {})
            pert_id = str(mm.get("pert_id", ""))
            pert_iname = str(mm.get("pert_iname", ""))
            cell_id = str(mm.get("cell_id", ""))
            cell_name = str(mm.get("cell_name", ""))
            moa = str(mm.get("moa", ""))
            target = str(mm.get("target", ""))

            e_series_all = expr_df.iloc[:, j]
            z_series = sp_by_tissue_dict[t]
            genes = z_series.index.intersection(e_series_all.index)
            if len(genes) < min_genes:
                continue

            z_sub = zscore_series(winsorize(z_series.loc[genes], q=winsor_q))
            e_sub = e_series_all.loc[genes]

            rho_all, p_s_all = spearmanr(z_sub.values, (-e_sub).values)
            r_all, p_p_all = pearsonr(z_sub.values, (-e_sub).values)
            method_p = [p_s_all, p_p_all]
            method_sign = [-rho_all, -r_all]

            for K in K_list:
                up_set, dn_set = select_gene_sets_from_Z(z_sub, K)
                subset_genes = list(set(list(up_set) + list(dn_set)))
                zs = z_sub.loc[subset_genes]
                es = e_sub.loc[subset_genes]

                rho_k, p_s_k = spearmanr(zs.values, (-es).values)
                r_k, p_p_k = pearsonr(zs.values, (-es).values)
                tau_k, p_k, cs_k, n_g = ks_connectivity_fast(z_sub, e_sub, K, alpha=alpha, nperm=nperm)

                method_p.extend([p_s_k, p_p_k, p_k])
                method_sign.extend([-rho_k, -r_k, -tau_k])

            Z_meta, p_meta = stouffer(method_p, signs=np.array(method_sign))
            final_rows.append({
                "tissue": t, "sig_id": sig_id, "pert_id": pert_id, "pert_iname": pert_iname,
                "cell_id": cell_id, "cell_name": cell_name, "moa": moa, "target": target,
                "Z_meta": Z_meta, "p_meta": p_meta
            })

            calc_cnt += 1
            if (calc_cnt % 500) == 0:
                print(f"[{ts()}] [DEBUG] {t}: refinement progress {calc_cnt}/{len(sel_sig_ids)}, elapsed={(time.time() - t_start) / 60:.2f} min", flush=True)

    # In fast mode, if top_pct < 1.0, append remaining Stage-1 records (here top_pct=1.0 -> non is empty)
    if mode == "fast" and len(non) > 0:
        for _, rec in non.iterrows():
            final_rows.append({
                "tissue": t, "sig_id": rec["sig_id"], "pert_id": rec.get("pert_id", ""),
                "pert_iname": rec.get("pert_iname", ""), "cell_id": rec.get("cell_id", ""),
                "cell_name": rec.get("cell_name", ""), "moa": rec.get("moa", ""),
                "target": rec.get("target", ""), "Z_meta": rec["Z_screen"], "p_meta": rec["p_screen"]
            })

    return t, final_rows, logmsg

# ================= CLI arguments =================
def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sp_dir", required=True, help="Directory of S-PrediXcan per-tissue CSV files")
    ap.add_argument("--gwas_id", required=True, help="GWAS identifier used in filenames and output folders")
    ap.add_argument("--sig_info_qc", required=True, help="Signature metadata table (QC-ed)")
    ap.add_argument("--cid_list", default="ALL", help="List of signature IDs to include (or ALL)")
    ap.add_argument("--gctx", required=True, help="Path to GCTX HDF5 matrix of L1000 signatures")
    ap.add_argument("--gene_info", required=True, help="Gene annotation table with pr_gene_id and pr_gene_symbol")
    ap.add_argument("--gene_info_landmark", default="", help="Landmark gene table if using landmark-only mode")
    ap.add_argument("--out_root", default=r"D:\Drug\HCC\output", help="Output root directory")
    ap.add_argument("--batch_size", type=int, default=1000, help="Batch size for signature processing")
    ap.add_argument("--workers", type=int, default=4, help="Number of parallel processes for Stage 2")
    ap.add_argument("--window_size", type=int, default=200, help="Column window size for HDF5 reads")
    ap.add_argument("--use_landmark_only", type=int, default=0, help="1 to restrict to landmark genes")
    ap.add_argument("--winsor_q", type=float, default=0.01, help="Winsorization quantile (two-sided)")
    ap.add_argument("--min_genes", type=int, default=100, help="Minimum overlapping genes required")
    ap.add_argument("--k_full", default="50,100,250,500", help="K list for core tissues")
    ap.add_argument("--nperm_liver", type=int, default=200, help="Permutations for Liver")
    ap.add_argument("--nperm_key_mid", type=int, default=200, help="Permutations for other key tissues")
    ap.add_argument("--k_fast", default="50,100", help="K list for non-core tissues")
    ap.add_argument("--nperm_fast", type=int, default=200, help="Permutations for fast mode")
    ap.add_argument("--top_pct", type=float, default=1.0, help="Top fraction to pass to Stage 2 (1.0 = all)")
    ap.add_argument("--alpha", type=float, default=1.0, help="KS weighting exponent (|e|^alpha)")
    ap.add_argument("--top_n", type=int, default=15, help="Top-N per tissue in final Table 2")
    ap.add_argument("--q_threshold", type=float, default=0.05, help="Global q threshold for optional downstream filters")
    ap.add_argument("--force_full", type=int, default=0, help="Force full mode for all tissues (1=yes)")
    return ap.parse_args()

# ================= Main pipeline =================
def main():
    print(f"[{ts()}] Script started...", flush=True)
    args = parse_args()
    out_dir = Path(args.out_root) / args.gwas_id
    out_dir.mkdir(parents=True, exist_ok=True)
    t_global = time.time()

    K_FAST = [int(x) for x in args.k_fast.split(",")]
    K_FULL = [int(x) for x in args.k_full.split(",")]
    TOP_PCT = float(args.top_pct)  # should be 1.0 to pass all to Stage 2
    TOP_N = int(args.top_n)

    manager = mp.Manager()
    h5_lock = manager.Lock()

    # 1) Load S-PrediXcan disease signatures
    pattern = os.path.join(args.sp_dir, f"{args.gwas_id}__PM__*.csv")
    sp_files = sorted(glob(pattern))
    if len(sp_files) == 0:
        raise SystemExit("[ERROR] No matching S-PrediXcan result files found.")

    sp_list = []
    for fp in sp_files:
        df = pd.read_csv(fp)
        m = re.search(r"__PM__([^\.]+)\.csv$", os.path.basename(fp))
        tissue = m.group(1)
        df["tissue"] = tissue
        df = df[["gene_name", "zscore", "tissue"]].dropna().rename(columns={"gene_name": "gene_symbol", "zscore": "Z"})
        sp_list.append(df)

    sp_all = pd.concat(sp_list, axis=0)
    sp_by_tissue = {t: d.set_index("gene_symbol")["Z"] for t, d in sp_all.groupby("tissue")}
    tissue_names = sorted(sp_by_tissue.keys())
    print(f"[{ts()}] [INFO] Number of tissues loaded: {len(tissue_names)}", flush=True)
    print(f"[{ts()}] [INFO] Core (8) tissues: {sorted(list(KEY_TISSUES_CORE))}", flush=True)

    # 2) Pre-index metadata dictionaries for speed
    sig_info = _read_table_auto(args.sig_info_qc)
    for c in ["sig_id","pert_iname","pert_id","cell_id","cell_name","moa","target"]:
        if c in sig_info.columns:
            sig_info[c] = sig_info[c].astype(str).fillna("")
    sig_info_dict = sig_info.set_index("sig_id").to_dict(orient="index")

    row_ids, col_ids, rid2idx, cid2idx = build_gctx_index(args.gctx)
    selected_cids = list(cid2idx.keys()) if args.cid_list == "ALL" else pd.read_csv(args.cid_list, header=None)[0].astype(str).tolist()
    selected_cids = [cid for cid in selected_cids if cid in cid2idx]
    print(f"[{ts()}] [INFO] Available signature count: {len(selected_cids)}", flush=True)

    gene_info = pd.read_csv(args.gene_info, sep="\t", low_memory=False)
    if int(args.use_landmark_only) == 1:
        lm = pd.read_csv(args.gene_info_landmark, sep="\t", low_memory=False)
        landmark_ids = set(lm["pr_gene_id"].astype(str))
        rid_selected = [rid for rid in gene_info["pr_gene_id"].astype(str) if rid in landmark_ids]
    else:
        rid_selected = list(gene_info["pr_gene_id"].astype(str))

    row_idx_list = [rid2idx[rid] for rid in rid_selected if rid in rid2idx]
    row_idx_list = np.unique(np.array(row_idx_list, dtype=np.int64)); row_idx_list.sort()
    rid2symbol_map = dict(zip(gene_info["pr_gene_id"].astype(str), gene_info["pr_gene_symbol"]))

    # ================= Stage 1: Vectorized screening =================
    screen_rows = []
    total_batches = int(np.ceil(len(selected_cids) / args.batch_size))

    print(f"[{ts()}] [INFO] [SCREEN] Aligning disease signature matrices across tissues...", flush=True)
    # 1) De-duplicated gene index
    all_genes_in_expr = sorted(list(set(rid2symbol_map.values())))

    # 2) Collect all tissues into a dict, avoiding reindex pitfalls
    z_dict = {}
    for t in tissue_names:
        s = sp_by_tissue[t]
        if s.index.has_duplicates:
            s = s.groupby(s.index).mean()
        z_dict[t] = s.reindex(all_genes_in_expr)

    # 3) Build a single DataFrame with unified gene index
    z_matrix_df = pd.DataFrame(z_dict, index=all_genes_in_expr)
    z_matrix_df = z_matrix_df.apply(lambda s: zscore_series(winsorize(s, q=args.winsor_q)))

    batch_id = 0
    for i in range(0, len(selected_cids), args.batch_size):
        batch_id += 1
        cid_chunk = selected_cids[i:i + args.batch_size]
        t_batch_start = time.time()

        col_idx_list = np.array([cid2idx[cid] for cid in cid_chunk], dtype=np.int64)
        expr_df = load_gctx_cols_windowed(args.gctx, row_idx_list, col_idx_list, row_ids, rid2symbol_map, window_size=args.window_size, lock=None)
        expr_df = expr_df.apply(lambda s: zscore_series(winsorize(s, q=args.winsor_q)))

        common_genes = z_matrix_df.index.intersection(expr_df.index)
        if len(common_genes) < args.min_genes:
            continue

        # Core acceleration: dot product across 49 tissues × N drugs
        Z_sub_mat = z_matrix_df.loc[common_genes].values
        E_sub_mat = expr_df.loc[common_genes].values

        n_g = len(common_genes)
        r_matrix = np.dot(Z_sub_mat.T, -E_sub_mat) / (n_g - 1)
        r_matrix = np.clip(r_matrix, -0.9999, 0.9999)
        t_stat = r_matrix * np.sqrt((n_g - 2) / (1 - r_matrix**2))
        p_matrix = 2 * norm.sf(np.abs(t_stat))

        # Pack results
        for s_idx, sig_id in enumerate(cid_chunk):
            mm = sig_info_dict.get(sig_id, {})
            for t_idx, t in enumerate(tissue_names):
                r_val = r_matrix[t_idx, s_idx]
                p_val = p_matrix[t_idx, s_idx]
                Z_screen, p_screen = stouffer([p_val, p_val], signs=[r_val, r_val])

                screen_rows.append({
                    "tissue": t, "sig_id": sig_id, "pert_id": mm.get("pert_id", ""),
                    "pert_iname": mm.get("pert_iname", ""), "cell_id": mm.get("cell_id", ""),
                    "cell_name": mm.get("cell_name", ""), "moa": mm.get("moa", ""),
                    "target": mm.get("target", ""), "Z_screen": Z_screen, "p_screen": p_screen
                })

        dt = time.time() - t_batch_start
        rate = len(cid_chunk) / dt if dt > 0 else 1
        remain_sig = len(selected_cids) - min(batch_id * args.batch_size, len(selected_cids))
        print(f"[{ts()}] [SCREEN] Batch {batch_id}/{total_batches} done | speed: {rate:.1f} sig/s | ETA Stage 1: {(remain_sig/rate)/60:.1f} min", flush=True)

    screen_df = pd.DataFrame(screen_rows)
    screen_df.to_csv(out_dir / "screen_pairs_all.csv", index=False)
    print(f"[{ts()}] [INFO] Stage 1 screening finished; results written.", flush=True)

    # ================= Stage 2: Multiprocess refinement per tissue =================
    tissue_args = []
    for t in tissue_names:
        df_t = screen_df[screen_df["tissue"] == t].copy()
        tissue_args.append((
            t, df_t, args.gctx, row_idx_list, row_ids, rid2symbol_map,
            K_FAST, args.nperm_fast, TOP_PCT, K_FULL, args.nperm_liver, args.nperm_key_mid,
            args.alpha, args.batch_size, args.winsor_q, args.min_genes,
            sig_info_dict, selected_cids, cid2idx, args.window_size, sp_by_tissue,
            int(args.force_full), h5_lock
        ))

    print(f"[{ts()}] [INFO] Launching Stage 2 refinement with multiprocessing; workers: {args.workers} ...", flush=True)
    final_rows = []
    with mp.Pool(processes=args.workers) as pool:
        for i, res in enumerate(pool.imap_unordered(refine_one_tissue, tissue_args), 1):
            t_name, rows, logmsg = res
            print(f"[{ts()}] [REFINE] Progress: {i}/{len(tissue_args)} tissues completed -> {t_name} | {logmsg}", flush=True)
            final_rows.extend(rows)

    pair_df = pd.DataFrame(final_rows)
    pair_df.to_csv(out_dir / f"pairwise_meta_hcc_parallel_{args.gwas_id}.csv", index=False)

    # ================= Multi-level aggregation and outputs =================
    # Tissue-level drug integration
    drug_tissue = []
    for (t, pert_id, pert_iname), sub in pair_df.groupby(["tissue", "pert_id", "pert_iname"]):
        Z_tissue, p_tissue = stouffer(sub["p_meta"].values, signs=sub["Z_meta"].values)
        drug_tissue.append({
            "tissue": t, "pert_id": pert_id, "pert_iname": pert_iname,
            "n_contexts": len(sub), "Z_agg_tissue": Z_tissue, "p_agg_tissue": p_tissue
        })
    drug_tissue_df = pd.DataFrame(drug_tissue)

    drug_tissue_final = []
    for t, sub in drug_tissue_df.groupby("tissue"):
        sub = sub.copy()
        sub["q_tissue"] = multipletests(sub["p_agg_tissue"].values, method="fdr_bh")[1]
        sub = sub.sort_values(["q_tissue", "p_agg_tissue"]).reset_index(drop=True)
        sub["Rank_in_tissue"] = np.arange(1, len(sub) + 1)
        drug_tissue_final.append(sub)
    drug_tissue_final = pd.concat(drug_tissue_final, axis=0)
    drug_tissue_final.to_csv(out_dir / f"drug_tissue_pq_hcc_parallel_{args.gwas_id}.csv", index=False)

    # Generate Table 2
    rep_cell = pair_df.sort_values("p_meta").groupby(["tissue", "pert_id", "pert_iname"], as_index=False).first()[["tissue", "pert_id", "pert_iname", "cell_name"]]

    def make_brief(row):
        moa, target = str(row.get('moa', '')).strip(), str(row.get('target', '')).strip()
        return f"{moa}; {target}" if moa and target else (moa if moa else target)

    brief_map = pair_df[["pert_id", "moa", "target"]].drop_duplicates("pert_id").copy()
    brief_map["brief"] = brief_map.apply(make_brief, axis=1)

    table2 = drug_tissue_final.merge(rep_cell, on=["tissue", "pert_id", "pert_iname"], how="left").merge(brief_map[["pert_id", "brief"]], on="pert_id", how="left")

    tb_list = []
    for t, sub in table2.groupby("tissue"):
        top_t = sub.sort_values(["q_tissue", "p_agg_tissue", "Rank_in_tissue"]).head(TOP_N).copy()
        top_t.rename(columns={
            "pert_iname": "Drug", "cell_name": "Cell line", "tissue": "Tissue",
            "Rank_in_tissue": "Rank in tissue", "p_agg_tissue": "P value",
            "q_tissue": "q value", "brief": "Brief description"
        }, inplace=True)
        tb_list.append(top_t[["Drug", "Cell line", "Tissue", "Rank in tissue", "P value", "q value", "Brief description", "n_contexts", "Z_agg_tissue"]])

    if tb_list:
        pd.concat(tb_list, axis=0).to_csv(out_dir / f"Table2_selected_repositioning_hits_49tissues_HCC_parallel_{args.gwas_id}.csv", index=False)

    # Global drug-level aggregation
    drug_level = []
    for pert_id, sub in drug_tissue_final.groupby("pert_id"):
        Z_drug, p_drug = stouffer(sub["p_agg_tissue"].values, signs=sub["Z_agg_tissue"].values)
        drug_level.append({
            "pert_id": pert_id, "pert_iname": sub["pert_iname"].iloc[0],
            "Z_agg_drug": Z_drug, "p_agg_drug": p_drug, "n_tissues": len(sub)
        })
    drug_level_df = pd.DataFrame(drug_level)
    drug_level_df["q_drug"] = multipletests(drug_level_df["p_agg_drug"].values, method="fdr_bh")[1]
    drug_level_df.sort_values("q_drug").to_csv(out_dir / f"drug_level_overall_HCC_parallel_{args.gwas_id}.csv", index=False)

    print(f"[{ts()}] === Pipeline finished successfully! Total time: {(time.time() - t_global) / 60:.2f} min ===", flush=True)

if __name__ == "__main__":
    try:
        mp.set_start_method("spawn", force=True)
    except RuntimeError:
        pass
    mp.freeze_support()
    main()