#!/usr/bin/env python3

import logging
import math
import re
import warnings
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
from matplotlib.patches import Patch
from scipy import sparse
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial import cKDTree
from scipy.spatial.distance import pdist
from scipy.stats import pearsonr


# Relative project paths
PROJECT_DIR = Path(".")
DATA_DIR = PROJECT_DIR / "data"
GSMAP_DIR = PROJECT_DIR / "gsmap_results"
OUTPUT_DIR = PROJECT_DIR / "figures"
COMBINED_CAUCHY_DIR = GSMAP_DIR / "cauchy_combined"

SAMPLES = {
    "slice1_mouse_brain": {
        "h5ad": DATA_DIR / "slice1.h5ad",
        "gsmap_dir": GSMAP_DIR / "slice1_mouse_brain",
    },
    "slice2_mouse_brain": {
        "h5ad": DATA_DIR / "slice2.h5ad",
        "gsmap_dir": GSMAP_DIR / "slice2_mouse_brain",
    },
}

TRAIT_KEYS = [
    "1_Mood_swings",
    "2_Miserableness",
    "3_Irritability",
    "4_Sensitivity",
    "5_Fed_up_feelings",
    "6_Nervous_feelings",
    "7_Anxious_feelings",
    "8_Tense",
    "9_Worry_too_long_after_embarrassment",
    "10_Suffer_from_nerves",
    "11_Loneliness",
    "12_Guilty_feelings",
    "13_Neuroticism_score_rint",
]
TRAIT_KEYS_NO_TOTAL = [
    trait for trait in TRAIT_KEYS if trait != "13_Neuroticism_score_rint"
]

TRAIT_LABELS = {
    "1_Mood_swings": "Mood swings",
    "2_Miserableness": "Miserableness",
    "3_Irritability": "Irritability",
    "4_Sensitivity": "Sensitivity",
    "5_Fed_up_feelings": "Fed-up feelings",
    "6_Nervous_feelings": "Nervous feelings",
    "7_Anxious_feelings": "Anxious feelings",
    "8_Tense": "Tense",
    "9_Worry_too_long_after_embarrassment": "Worry after embarrassment",
    "10_Suffer_from_nerves": "Suffer from nerves",
    "11_Loneliness": "Loneliness",
    "12_Guilty_feelings": "Guilty feelings",
    "13_Neuroticism_score_rint": "Neuroticism score",
}

MAJOR_ORDER = ["GLU", "GABA", "OPC", "Oligo", "Micro", "Astro", "Epen", "Endo", "CHOR"]
CELLTYPE_GROUPS = {
    "GLU": [
        "DIME_N_GLU", "L2/3_IT_GLU", "L4/5_IT_GLU", "L5_IT_GLU", "L5_PT_GLU",
        "L5/6_NP_GLU", "L6_N_GLU", "L6_IT_GLU", "L6_Car3_GLU", "TE_N_GLU",
        "DG_N_GLU", "CA1_N_GLU", "CA2_N_GLU", "CA3_N_GLU", "RH_N_GLU",
    ],
    "GABA": [
        "CNU_N_GABA", "DIME_N_GABA", "TE_N_GABA", "TE_N_GABA_SST",
        "TE_N_GABA_PVALB", "TE_N_GABA_LAMP5", "TE_N_GABA_SNCG",
        "TE_N_GABA_VIP", "TE_N_GABA_RELN", "TE_N_GABA_CHODL",
        "OB_N_GABA", "RH_N_GABA",
    ],
    "OPC": ["OPC"],
    "Oligo": ["OL"],
    "Micro": ["MGL"],
    "Astro": ["ASC"],
    "Epen": ["EPC"],
    "Endo": ["EDC"],
    "CHOR": ["CHOR"],
}
FINE_TO_MAJOR = {
    fine: major for major, fine_types in CELLTYPE_GROUPS.items() for fine in fine_types
}

REGION_SOURCE = "combined_prefer"
MANTEL_SAMPLE_N = 300
MANTEL_PERMUTATIONS = 999
MANTEL_RANDOM_STATE = 20260310
K_SPATIAL = 8
FIG_DPI = 400
LOGP_CLIP_MIN = 1e-300
CENTER_LOGP = 5.0
STAR_CUTS = [(1e-15, "***"), (1e-11, "**"), (1e-7, "*")]

logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
warnings.filterwarnings("ignore")

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "DejaVu Sans"
mpl.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial", "Liberation Sans"]
mpl.rcParams["mathtext.fontset"] = "dejavusans"
mpl.rcParams["axes.linewidth"] = 0.8
mpl.rcParams["xtick.major.width"] = 0.8
mpl.rcParams["ytick.major.width"] = 0.8

cmap_pearson = LinearSegmentedColormap.from_list(
    "pearson_map", ["#35978f", "#f5f5f5", "#bf812d"], N=256
)
cmap_mantel = LinearSegmentedColormap.from_list(
    "mantel_map", ["#c7eae5", "#80cdc1", "#35978f"], N=256
)
cmap_bwr = LinearSegmentedColormap.from_list(
    "bwr_custom", ["#2b8cbe", "#ffffff", "#b30000"], N=256
)
cmap_bwr.set_bad("#f2f2f2")


def sanitize_key(s):
    s = re.sub(r"[^0-9A-Za-z_]+", "_", str(s))
    return re.sub(r"_+", "_", s).strip("_")

def safe_neglog10(x):
    x = pd.to_numeric(x, errors="coerce").astype(float)
    y = -np.log10(np.clip(x, LOGP_CLIP_MIN, 1.0))
    y[pd.isna(x)] = np.nan
    return y

def p_to_star(p):
    if pd.isna(p):
        return ""
    for cutoff, star in STAR_CUTS:
        if p < cutoff:
            return star
    return ""

def save_df(df, path):
    df.to_csv(path, sep="\t", index=True)

def get_trait_obs_cols(trait):
    return f"gsmap_p__{sanitize_key(trait)}", f"gsmap_nlog10p__{sanitize_key(trait)}"

def ensure_trait_series(adata, trait):
    p_col, nlog_col = get_trait_obs_cols(trait)
    p = pd.to_numeric(adata.obs[p_col], errors="coerce").astype(float)
    nlog = pd.to_numeric(adata.obs[nlog_col], errors="coerce").astype(float) if nlog_col in adata.obs.columns else safe_neglog10(p)
    return p, nlog

def get_xy(adata):
    if {"x_rot_um", "y_rot_um"}.issubset(adata.obs.columns):
        return (
            pd.to_numeric(adata.obs["x_rot_um"], errors="coerce").to_numpy(),
            pd.to_numeric(adata.obs["y_rot_um"], errors="coerce").to_numpy(),
        )
    if "spatial" in adata.obsm:
        arr = np.asarray(adata.obsm["spatial"])
        return arr[:, 0], arr[:, 1]
    return (
        pd.to_numeric(adata.obs["x"], errors="coerce").to_numpy(),
        pd.to_numeric(adata.obs["y"], errors="coerce").to_numpy(),
    )

def infer_point_size(n):
    return max(0.12, min(0.7, 35000 / max(n, 1)))

def dictlike_uns_to_df(obj):
    if isinstance(obj, pd.DataFrame):
        return obj.copy()
    cols = list(obj["columns"])
    return pd.DataFrame({c: np.array(obj["data"][c]) for c in cols})

def corr_fast(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a - a.mean()
    b = b - b.mean()
    den = np.sqrt(np.dot(a, a) * np.dot(b, b))
    return np.nan if den == 0 else float(np.dot(a, b) / den)

def safe_pearsonr(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 3 or np.nanstd(x) == 0 or np.nanstd(y) == 0:
        return np.nan, np.nan
    return pearsonr(x, y)

def mantel_test_1d(x, y, permutations=999, random_state=0):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 5:
        return np.nan, np.nan
    dx = pdist(x[:, None], metric="euclidean")
    dy = pdist(y[:, None], metric="euclidean")
    obs = corr_fast(dx, dy)
    rng = np.random.default_rng(random_state)
    cnt = 1
    for _ in range(permutations):
        yp = rng.permutation(y)
        rp = corr_fast(dx, pdist(yp[:, None], metric="euclidean"))
        if np.abs(rp) >= np.abs(obs):
            cnt += 1
    return obs, cnt / (permutations + 1)

def proportional_sample(df, n_total=300, group_col="sample", random_state=0):
    rng = np.random.default_rng(random_state)
    vc = df[group_col].value_counts()
    props = vc / vc.sum()
    alloc = (props * n_total).round().astype(int)
    diff = n_total - alloc.sum()
    order = props.sort_values(ascending=False).index.tolist()
    i = 0
    while diff != 0:
        g = order[i % len(order)]
        alloc[g] += 1 if diff > 0 else -1
        diff += -1 if diff > 0 else 1
        i += 1
    out = []
    for g, n in alloc.items():
        sub = df[df[group_col] == g]
        n = min(n, len(sub))
        idx = rng.choice(sub.index.to_numpy(), size=n, replace=False)
        out.append(sub.loc[idx])
    return pd.concat(out).sample(frac=1.0, random_state=random_state)

def make_diverging_norm(values, center=CENTER_LOGP):
    arr = np.asarray(values, float)
    arr = arr[np.isfinite(arr)]
    if len(arr) == 0:
        return Normalize(vmin=0, vmax=1)
    vmax = max(float(np.nanmax(arr)), 1.0)
    vcenter = center if vmax > center else max(vmax * 0.5, 1e-6)
    return TwoSlopeNorm(vmin=0.0, vcenter=vcenter, vmax=vmax)

def acat_cauchy(pvals):
    p = pd.to_numeric(pd.Series(pvals), errors="coerce").dropna().astype(float).values
    p = p[(p > 0) & (p <= 1)]
    if len(p) == 0:
        return np.nan
    p = np.clip(p, 1e-300, 1 - 1e-16)
    t = np.tan((0.5 - p) * np.pi)
    T = np.mean(t)
    if T > 1e15:
        out = 1.0 / (np.pi * T)
    elif T < -1e15:
        out = 1.0
    else:
        out = 0.5 - np.arctan(T) / np.pi
    return float(np.clip(out, 1e-300, 1.0))

def build_cell_trait_matrix(adata1, adata2, trait_keys):
    blocks = []
    for sample, adata in [("slice1_mouse_brain", adata1), ("slice2_mouse_brain", adata2)]:
        df = pd.DataFrame(index=adata.obs_names.astype(str))
        df["sample"] = sample
        for trait in trait_keys:
            _, nlog = ensure_trait_series(adata, trait)
            df[TRAIT_LABELS[trait]] = nlog.values
        blocks.append(df)
    out = pd.concat(blocks, axis=0)
    cols = [TRAIT_LABELS[t] for t in trait_keys]
    return out.dropna(subset=cols, how="any").copy()

def compute_pearson_mats(df):
    cols = list(df.columns)
    r = pd.DataFrame(np.eye(len(cols)), index=cols, columns=cols, dtype=float)
    p = pd.DataFrame(np.zeros((len(cols), len(cols))), index=cols, columns=cols, dtype=float)
    for i, c1 in enumerate(cols):
        for j, c2 in enumerate(cols):
            if j < i:
                continue
            if i != j:
                rv, pv = safe_pearsonr(df[c1].values, df[c2].values)
                r.loc[c1, c2] = r.loc[c2, c1] = rv
                p.loc[c1, c2] = p.loc[c2, c1] = pv
    return r, p

def compute_mantel_mats(df, sample_n=300, permutations=999, random_state=0):
    sampled = proportional_sample(df.copy(), n_total=sample_n, group_col="sample", random_state=random_state)
    X = sampled[[c for c in sampled.columns if c != "sample"]].copy()
    cols = list(X.columns)
    r = pd.DataFrame(np.eye(len(cols)), index=cols, columns=cols, dtype=float)
    p = pd.DataFrame(np.zeros((len(cols), len(cols))), index=cols, columns=cols, dtype=float)
    for i, c1 in enumerate(cols):
        for j, c2 in enumerate(cols):
            if j < i:
                continue
            if i != j:
                rv, pv = mantel_test_1d(X[c1].values, X[c2].values, permutations=permutations, random_state=random_state + i * 1000 + j)
                r.loc[c1, c2] = r.loc[c2, c1] = rv
                p.loc[c1, c2] = p.loc[c2, c1] = pv
    return r, p, sampled


def load_gsmap_results(sample_name, config, trait_keys):
    """Load an AnnData object and attach core GsMap post-analysis results.

    Spot-level files are matched to ``adata.obs_names`` through the first file
    column. Region-level Cauchy results are stored in ``adata.uns`` for later
    heatmap construction. Existing columns or ``uns`` entries are retained when
    the corresponding result file is not present.
    """
    adata = sc.read_h5ad(config["h5ad"])
    result_dir = Path(config["gsmap_dir"])

    for trait in trait_keys:
        trait_key = sanitize_key(trait)
        p_col, nlog_col = get_trait_obs_cols(trait)

        spot_file = (
            result_dir
            / "spatial_ldsc"
            / f"{sample_name}_{trait}.csv.gz"
        )
        if spot_file.exists():
            spot = pd.read_csv(spot_file, compression="infer").copy()
            spot = spot.rename(columns={spot.columns[0]: "row_id"})
            if "p" not in spot.columns:
                raise KeyError(f"Missing 'p' column in {spot_file}")

            spot_p = pd.Series(
                pd.to_numeric(spot["p"], errors="coerce").to_numpy(),
                index=spot["row_id"].astype(str),
            )
            spot_p = spot_p[~spot_p.index.duplicated(keep="first")]
            adata.obs[p_col] = adata.obs_names.astype(str).map(spot_p)
        elif p_col not in adata.obs.columns:
            raise FileNotFoundError(
                f"Neither {spot_file} nor adata.obs['{p_col}'] is available."
            )

        if nlog_col not in adata.obs.columns:
            adata.obs[nlog_col] = safe_neglog10(adata.obs[p_col]).to_numpy()

        cauchy_file = (
            result_dir
            / "cauchy_combination"
            / f"{sample_name}_{trait}.Cauchy.csv.gz"
        )
        uns_key = f"gsmap_cauchy__{trait_key}"
        if cauchy_file.exists():
            cauchy = pd.read_csv(cauchy_file, compression="infer").copy()
            if "annotation" not in cauchy.columns or "p_cauchy" not in cauchy.columns:
                raise KeyError(
                    f"Required columns 'annotation' and 'p_cauchy' are missing in {cauchy_file}"
                )
            cauchy["annotation"] = cauchy["annotation"].astype(str)
            cauchy["p_cauchy"] = pd.to_numeric(cauchy["p_cauchy"], errors="coerce")
            adata.uns[uns_key] = cauchy[["annotation", "p_cauchy"]].copy()

            region_p = dict(zip(cauchy["annotation"], cauchy["p_cauchy"]))
            region_p_col = f"gsmap_region_p_cauchy__{trait_key}"
            region_nlog_col = f"gsmap_region_nlog10p_cauchy__{trait_key}"
            adata.obs[region_p_col] = adata.obs["region"].astype(str).map(region_p)
            adata.obs[region_nlog_col] = safe_neglog10(adata.obs[region_p_col]).to_numpy()
        elif uns_key not in adata.uns:
            warnings.warn(f"Region-level Cauchy result not found for {sample_name}: {trait}")

    return adata


def read_combined_cauchy_file(trait, combined_dir=COMBINED_CAUCHY_DIR):
    """Read an optional Cauchy result combined across spatial samples."""
    path = Path(combined_dir) / f"{trait}.combined.cauchy.csv.gz"
    if path.exists():
        return pd.read_csv(path, compression="infer"), path
    return None, None


def get_region_table(
    adata_ref,
    trait,
    source="combined_prefer",
    combined_dir=COMBINED_CAUCHY_DIR,
):
    """Retrieve region-level GsMap Cauchy results from AnnData or a combined file."""
    combined_key = f"gsmap_combined_cauchy__{sanitize_key(trait)}"
    single_key = f"gsmap_cauchy__{sanitize_key(trait)}"

    if source == "combined_only":
        if combined_key in adata_ref.uns:
            return dictlike_uns_to_df(adata_ref.uns[combined_key]), "uns_combined"
        table, _ = read_combined_cauchy_file(trait, combined_dir)
        if table is not None:
            return table, "file_combined"
        raise KeyError(f"No combined region-level result found for {trait}")

    if source == "single_only":
        if single_key not in adata_ref.uns:
            raise KeyError(f"No sample-level region result found for {trait}")
        return dictlike_uns_to_df(adata_ref.uns[single_key]), "uns_single"

    if combined_key in adata_ref.uns:
        return dictlike_uns_to_df(adata_ref.uns[combined_key]), "uns_combined"

    table, _ = read_combined_cauchy_file(trait, combined_dir)
    if table is not None:
        return table, "file_combined"

    if single_key in adata_ref.uns:
        return dictlike_uns_to_df(adata_ref.uns[single_key]), "uns_single"

    raise KeyError(f"No region-level result found for {trait}")


def cluster_rows(df):
    x = df.copy().fillna(0.0)
    if x.shape[0] <= 2:
        return df
    order = leaves_list(linkage(x.values, method="average", metric="euclidean"))
    return df.iloc[order, :]

def build_region_trait_mats(adata_ref, trait_keys, source="combined_prefer", combined_dir=COMBINED_CAUCHY_DIR):
    p_list, source_rows = [], []
    for trait in trait_keys:
        df, src = get_region_table(adata_ref, trait, source=source, combined_dir=combined_dir)
        df = df.copy()
        df["annotation"] = df["annotation"].astype(str)
        df["p_cauchy"] = pd.to_numeric(df["p_cauchy"], errors="coerce")
        p_list.append(df.set_index("annotation")["p_cauchy"].rename(TRAIT_LABELS[trait]))
        source_rows.append({"trait": trait, "source": src})
    p_mat = pd.concat(p_list, axis=1).drop(index="root", errors="ignore")
    nlog_mat = p_mat.apply(safe_neglog10)
    nlog_mat = cluster_rows(nlog_mat)
    p_mat = p_mat.loc[nlog_mat.index, nlog_mat.columns]
    return p_mat, nlog_mat, pd.DataFrame(source_rows)

def get_major_series(adata):
    s = adata.obs["fine"].astype(str).map(FINE_TO_MAJOR)
    return pd.Series(pd.Categorical(s, categories=MAJOR_ORDER, ordered=True), index=adata.obs_names.astype(str), name="major")

def build_trait_major_mats(adata1, adata2, trait_keys):
    blocks = []
    for sample, adata in [("slice1_mouse_brain", adata1), ("slice2_mouse_brain", adata2)]:
        df = pd.DataFrame(index=adata.obs_names.astype(str))
        df["sample"] = sample
        df["major"] = get_major_series(adata).astype(object)
        for trait in trait_keys:
            p, _ = ensure_trait_series(adata, trait)
            df[trait] = p.values
        blocks.append(df)
    all_df = pd.concat(blocks, axis=0)
    all_df = all_df[all_df["major"].notna()].copy()
    counts = all_df["major"].value_counts().reindex(MAJOR_ORDER).dropna().astype(int)

    p_rows = []
    for trait in trait_keys:
        row = {major: acat_cauchy(all_df.loc[all_df["major"] == major, trait].values) for major in MAJOR_ORDER}
        p_rows.append(pd.Series(row, name=TRAIT_LABELS[trait]))
    p_mat = pd.DataFrame(p_rows).reindex(index=[TRAIT_LABELS[t] for t in trait_keys], columns=MAJOR_ORDER)
    nlog_mat = p_mat.apply(safe_neglog10)
    return p_mat, nlog_mat, counts

def build_knn_w(coords, k=8):
    n = coords.shape[0]
    k_eff = min(k, max(n - 1, 1))
    tree = cKDTree(coords)
    _, idx = tree.query(coords, k=k_eff + 1)
    rows = np.repeat(np.arange(n), k_eff)
    cols = idx[:, 1:k_eff + 1].reshape(-1)
    data = np.ones(len(rows), dtype=float)
    W = sparse.csr_matrix((data, (rows, cols)), shape=(n, n))
    W = ((W + W.T) > 0).astype(float).tocsr()
    rs = np.asarray(W.sum(axis=1)).ravel()
    inv = np.divide(1.0, rs, out=np.zeros_like(rs), where=rs > 0)
    W = sparse.diags(inv) @ W
    return W.tocsr()

def zscore_cols(X, ddof=0):
    X = np.asarray(X, float)
    mu = np.nanmean(X, axis=0)
    sd = np.nanstd(X, axis=0, ddof=ddof)
    sd[sd == 0] = 1.0
    return (X - mu) / sd

def get_spatial_trait_block(adata, trait_keys):
    x, y = get_xy(adata)
    vals = []
    for trait in trait_keys:
        _, nlog = ensure_trait_series(adata, trait)
        vals.append(pd.to_numeric(nlog, errors="coerce").to_numpy())
    X = np.column_stack(vals)
    coords = np.column_stack([x, y])
    mask = np.isfinite(coords).all(axis=1) & np.isfinite(X).all(axis=1)
    return X[mask], coords[mask]

def build_spatial_trait_mats(adata1, adata2, trait_keys, k=8):
    X1, C1 = get_spatial_trait_block(adata1, trait_keys)
    X2, C2 = get_spatial_trait_block(adata2, trait_keys)
    W1 = build_knn_w(C1, k=k)
    W2 = build_knn_w(C2, k=k)
    W = sparse.block_diag([W1, W2], format="csr")
    X = np.vstack([X1, X2])

    Z_lee = zscore_cols(X, ddof=0)
    WZ = W @ Z_lee
    lee = (WZ.T @ WZ) / W.shape[0]
    lee = pd.DataFrame(lee, index=[TRAIT_LABELS[t] for t in trait_keys], columns=[TRAIT_LABELS[t] for t in trait_keys])

    Z_moran = zscore_cols(X, ddof=1)
    moran = (Z_moran.T @ (W @ Z_moran)) / (W.shape[0] - 1.0)
    moran = pd.DataFrame(moran, index=[TRAIT_LABELS[t] for t in trait_keys], columns=[TRAIT_LABELS[t] for t in trait_keys])

    np.fill_diagonal(lee.values, np.nan)
    np.fill_diagonal(moran.values, np.nan)
    return lee, moran

def draw_heatmap(ax, value_df, p_df=None, cmap=None, norm=None, vmin=None, vmax=None, xrot=45, show_stars=True, ylabels=True, xfontsize=10, yfontsize=10):
    arr = np.ma.masked_invalid(value_df.to_numpy(dtype=float))
    im = ax.imshow(arr, cmap=cmap, norm=norm, vmin=vmin, vmax=vmax, aspect="auto")

    nr, nc = value_df.shape
    ax.set_xticks(np.arange(nc))
    ax.set_yticks(np.arange(nr))
    ax.set_xticklabels(list(value_df.columns), rotation=xrot, ha="right" if xrot else "center", fontsize=xfontsize)
    ax.set_yticklabels(list(value_df.index) if ylabels else [""] * nr, fontsize=yfontsize)

    ax.set_xticks(np.arange(-0.5, nc, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, nr, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.0)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.tick_params(length=0)

    for spine in ax.spines.values():
        spine.set_visible(False)

    if show_stars and p_df is not None:
        for i in range(nr):
            for j in range(nc):
                s = p_to_star(p_df.iloc[i, j])
                if s:
                    ax.text(j, i, s, ha="center", va="center", fontsize=8, color="black")
    return im

def plot_spatial(adata, trait, out_png, out_pdf, cmap):
    x, y = get_xy(adata)
    _, nlog = ensure_trait_series(adata, trait)
    v = nlog.to_numpy(dtype=float)
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(v)

    fig, ax = plt.subplots(figsize=(6.0, 6.0))
    sca = ax.scatter(
        x[m], y[m], c=v[m], s=infer_point_size(len(adata)),
        cmap=cmap, norm=make_diverging_norm(v[m]), linewidths=0, rasterized=True
    )
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(sca, ax=ax, fraction=0.046, pad=0.02)
    cbar.set_label("-log10(p)", rotation=270, labelpad=16)
    cbar.outline.set_visible(False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=FIG_DPI, bbox_inches="tight")
    fig.savefig(out_pdf, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)

def plot_spatial_overview(adata, trait_keys, out_png, out_pdf, cmap):
    n, ncols = len(trait_keys), 4
    nrows = math.ceil(n / ncols)
    x, y = get_xy(adata)

    fig, axes = plt.subplots(nrows, ncols, figsize=(4.0 * ncols, 4.0 * nrows))
    axes = np.array(axes).reshape(nrows, ncols)

    for idx, trait in enumerate(trait_keys):
        ax = axes.flat[idx]
        _, nlog = ensure_trait_series(adata, trait)
        v = nlog.to_numpy(dtype=float)
        m = np.isfinite(x) & np.isfinite(y) & np.isfinite(v)
        sca = ax.scatter(
            x[m], y[m], c=v[m], s=infer_point_size(len(adata)),
            cmap=cmap, norm=make_diverging_norm(v[m]), linewidths=0, rasterized=True
        )
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        cbar = fig.colorbar(sca, ax=ax, fraction=0.045, pad=0.01)
        cbar.set_label("-log10(p)", rotation=270, labelpad=12)
        cbar.outline.set_visible(False)

    for idx in range(n, nrows * ncols):
        axes.flat[idx].axis("off")

    fig.tight_layout()
    fig.savefig(out_png, dpi=FIG_DPI, bbox_inches="tight")
    fig.savefig(out_pdf, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)

def build_region_palette(adata1, adata2):
    palette = {}

    def add_from_adata(adata):
        if "region_colors" in adata.uns and hasattr(adata.obs["region"], "cat"):
            cats = list(adata.obs["region"].cat.categories)
            cols = list(adata.uns["region_colors"])
            if len(cats) == len(cols):
                for c, col in zip(cats, cols):
                    palette[str(c)] = col

    add_from_adata(adata1)
    add_from_adata(adata2)

    union = list(pd.Index(adata1.obs["region"].astype(str)).union(pd.Index(adata2.obs["region"].astype(str))))
    missing = [r for r in union if r not in palette]
    if missing:
        cm = plt.get_cmap("tab20", len(missing))
        for i, r in enumerate(missing):
            palette[r] = mpl.colors.to_hex(cm(i))
    return palette

def plot_region_pair(adata1, adata2, palette, out_png, out_pdf):
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    for ax, adata in zip(axes, [adata1, adata2]):
        x, y = get_xy(adata)
        reg = adata.obs["region"].astype(str)
        colors = reg.map(palette).fillna("#bdbdbd").to_numpy()
        m = np.isfinite(x) & np.isfinite(y)
        ax.scatter(x[m], y[m], c=colors[m], s=infer_point_size(len(adata)), linewidths=0, rasterized=True)
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=FIG_DPI, bbox_inches="tight")
    fig.savefig(out_pdf, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)

def plot_region_legend(palette, adata1, adata2, out_png, out_pdf):
    union = list(pd.Index(adata1.obs["region"].astype(str)).union(pd.Index(adata2.obs["region"].astype(str))))
    labels = [x for x in union if x != "root"]
    handles = [Patch(facecolor=palette[x], edgecolor="none", label=x) for x in labels]
    fig, ax = plt.subplots(figsize=(12, max(2.5, 0.32 * len(labels))))
    ax.axis("off")
    ax.legend(handles=handles, loc="center", ncol=4, frameon=False, fontsize=10, handlelength=1.2, columnspacing=1.2)
    fig.tight_layout()
    fig.savefig(out_png, dpi=FIG_DPI, bbox_inches="tight")
    fig.savefig(out_pdf, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)

def plot_dual_heatmap(left_df, right_df, right_p_df, left_label, right_label, left_cmap, right_cmap, out_png, out_pdf):
    fig = plt.figure(figsize=(15.2, 8.4))
    gs = fig.add_gridspec(1, 4, width_ratios=[0.05, 1.05, 0.86, 0.05], wspace=0.18)

    cax_left = fig.add_subplot(gs[0, 0])
    ax_left = fig.add_subplot(gs[0, 1])
    ax_right = fig.add_subplot(gs[0, 2], sharey=ax_left)
    cax_right = fig.add_subplot(gs[0, 3])

    im_left = draw_heatmap(
        ax_left,
        left_df,
        p_df=None,
        cmap=left_cmap,
        norm=Normalize(vmin=-1, vmax=1),
        xrot=45,
        show_stars=False,
        ylabels=True,
    )
    im_right = draw_heatmap(
        ax_right,
        right_df,
        p_df=right_p_df,
        cmap=right_cmap,
        norm=make_diverging_norm(right_df.to_numpy()),
        xrot=0,
        show_stars=True,
        ylabels=False,
    )

    cb_left = fig.colorbar(im_left, cax=cax_left)
    cb_left.set_label(left_label, rotation=90, labelpad=16)
    cb_left.outline.set_visible(False)

    cb_right = fig.colorbar(im_right, cax=cax_right)
    cb_right.set_label(right_label, rotation=270, labelpad=16)
    cb_right.outline.set_visible(False)

    ax_left.set_xlabel("")
    ax_left.set_ylabel("")
    ax_right.set_xlabel("")
    ax_right.set_ylabel("")
    plt.setp(ax_right.get_yticklabels(), visible=False)

    fig.tight_layout()
    fig.savefig(out_png, dpi=FIG_DPI, bbox_inches="tight")
    fig.savefig(out_pdf, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)

def plot_dual_heatmap_transposed(left_df, left_p_df, right_df, right_p_df, left_label, right_label, left_cmap, right_cmap, out_png, out_pdf):
    fig = plt.figure(figsize=(18.0, 8.4))
    gs = fig.add_gridspec(1, 4, width_ratios=[0.055, 1.28, 0.88, 0.045], wspace=0.24)

    cax_left = fig.add_subplot(gs[0, 0])
    ax_left = fig.add_subplot(gs[0, 1])
    ax_right = fig.add_subplot(gs[0, 2], sharey=ax_left)
    cax_right = fig.add_subplot(gs[0, 3])

    im_left = draw_heatmap(
        ax_left,
        left_df,
        p_df=left_p_df,
        cmap=left_cmap,
        norm=make_diverging_norm(left_df.to_numpy()),
        xrot=45,
        show_stars=True,
        ylabels=True,
        xfontsize=9,
        yfontsize=10,
    )
    im_right = draw_heatmap(
        ax_right,
        right_df,
        p_df=right_p_df,
        cmap=right_cmap,
        norm=make_diverging_norm(right_df.to_numpy()),
        xrot=0,
        show_stars=True,
        ylabels=False,
        xfontsize=10,
        yfontsize=10,
    )

    cb_left = fig.colorbar(im_left, cax=cax_left)
    cb_left.set_label(left_label, rotation=90, labelpad=20)
    cb_left.outline.set_visible(False)

    cb_right = fig.colorbar(im_right, cax=cax_right)
    cb_right.set_label(right_label, rotation=270, labelpad=16)
    cb_right.outline.set_visible(False)

    ax_left.set_xlabel("")
    ax_left.set_ylabel("")
    ax_right.set_xlabel("")
    ax_right.set_ylabel("")
    plt.setp(ax_right.get_yticklabels(), visible=False)

    fig.tight_layout()
    fig.savefig(out_png, dpi=FIG_DPI, bbox_inches="tight")
    fig.savefig(out_pdf, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


def create_output_directories():
    subdirectories = [
        "01_trait_correlation",
        "02_region_trait_heatmap",
        "03_spatial_maps/slice1_mouse_brain",
        "03_spatial_maps/slice2_mouse_brain",
        "04_major_celltype_heatmap",
        "05_slice2_HPF_spatial_maps",
        "06_lee_with_major_heatmap",
        "07_moran_bv_with_major_heatmap",
        "08_region_maps",
        "09_region_transposed_with_major_heatmap",
    ]
    for subdirectory in subdirectories:
        (OUTPUT_DIR / subdirectory).mkdir(parents=True, exist_ok=True)


def run_analysis():
    create_output_directories()

    adatas = {
        sample_name: load_gsmap_results(sample_name, config, TRAIT_KEYS)
        for sample_name, config in SAMPLES.items()
    }
    adata1 = adatas["slice1_mouse_brain"]
    adata2 = adatas["slice2_mouse_brain"]

    # 1. Trait-level Pearson and Mantel correlations.
    cell_trait_df = build_cell_trait_matrix(adata1, adata2, TRAIT_KEYS_NO_TOTAL)
    correlation_columns = [TRAIT_LABELS[trait] for trait in TRAIT_KEYS_NO_TOTAL]
    pearson_r, pearson_p = compute_pearson_mats(cell_trait_df[correlation_columns])
    mantel_r, mantel_p, _ = compute_mantel_mats(
        cell_trait_df[["sample"] + correlation_columns],
        sample_n=MANTEL_SAMPLE_N,
        permutations=MANTEL_PERMUTATIONS,
        random_state=MANTEL_RANDOM_STATE,
    )

    correlation_dir = OUTPUT_DIR / "01_trait_correlation"
    save_df(pearson_r, correlation_dir / "trait_pearson_r.tsv")
    save_df(pearson_p, correlation_dir / "trait_pearson_p.tsv")
    save_df(mantel_r, correlation_dir / "trait_mantel_r.tsv")
    save_df(mantel_p, correlation_dir / "trait_mantel_p.tsv")

    fig, axes = plt.subplots(1, 2, figsize=(17, 8))
    im1 = draw_heatmap(
        axes[0], pearson_r, pearson_p, cmap=cmap_pearson,
        vmin=-1, vmax=1, xrot=45, show_stars=True, ylabels=True,
    )
    im2 = draw_heatmap(
        axes[1], mantel_r, mantel_p, cmap=cmap_mantel,
        vmin=0, vmax=1, xrot=45, show_stars=True, ylabels=True,
    )
    cb1 = fig.colorbar(im1, ax=axes[0], fraction=0.046, pad=0.04)
    cb1.set_label("Pearson's r", rotation=270, labelpad=18)
    cb1.outline.set_visible(False)
    cb2 = fig.colorbar(im2, ax=axes[1], fraction=0.046, pad=0.04)
    cb2.set_label("Mantel's r", rotation=270, labelpad=18)
    cb2.outline.set_visible(False)
    fig.tight_layout()
    fig.savefig(correlation_dir / "trait_correlation_heatmaps.png", dpi=FIG_DPI, bbox_inches="tight")
    fig.savefig(correlation_dir / "trait_correlation_heatmaps.pdf", dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)

    # 2. Region-by-trait Cauchy-combination heatmap.
    region_p, region_nlog, region_source = build_region_trait_mats(
        adata1,
        TRAIT_KEYS,
        source=REGION_SOURCE,
        combined_dir=COMBINED_CAUCHY_DIR,
    )
    region_dir = OUTPUT_DIR / "02_region_trait_heatmap"
    save_df(region_p, region_dir / "region_trait_p_cauchy.tsv")
    save_df(region_nlog, region_dir / "region_trait_nlog10_p_cauchy.tsv")
    region_source.to_csv(region_dir / "region_trait_source.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(10.8, 9.0))
    im = draw_heatmap(
        ax, region_nlog, region_p, cmap=cmap_bwr,
        norm=make_diverging_norm(region_nlog.to_numpy()),
        xrot=45, show_stars=True, ylabels=True,
    )
    cb = fig.colorbar(im, ax=ax, fraction=0.032, pad=0.02)
    cb.set_label("-log10(p_cauchy)", rotation=270, labelpad=18)
    cb.outline.set_visible(False)
    fig.tight_layout()
    fig.savefig(region_dir / "region_trait_nlog10_p_cauchy_heatmap.png", dpi=FIG_DPI, bbox_inches="tight")
    fig.savefig(region_dir / "region_trait_nlog10_p_cauchy_heatmap.pdf", dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)

    # 3. Spot-level GsMap spatial maps for both slices.
    for sample_name, adata in adatas.items():
        sample_dir = OUTPUT_DIR / "03_spatial_maps" / sample_name
        for trait in TRAIT_KEYS:
            tag = sanitize_key(trait)
            plot_spatial(
                adata,
                trait,
                sample_dir / f"{sample_name}.{tag}.png",
                sample_dir / f"{sample_name}.{tag}.pdf",
                cmap_bwr,
            )

    # 4. Major-cell-type Cauchy-combination heatmap.
    trait_major_p, trait_major_nlog, major_counts = build_trait_major_mats(
        adata1, adata2, TRAIT_KEYS_NO_TOTAL
    )
    major_dir = OUTPUT_DIR / "04_major_celltype_heatmap"
    save_df(trait_major_p, major_dir / "trait_major_p_cauchy.tsv")
    save_df(trait_major_nlog, major_dir / "trait_major_nlog10_p_cauchy.tsv")
    major_counts.to_csv(major_dir / "major_celltype_counts.tsv", sep="\t", header=True)

    fig, ax = plt.subplots(figsize=(9.5, 8.2))
    im = draw_heatmap(
        ax, trait_major_nlog, trait_major_p, cmap=cmap_pearson,
        norm=make_diverging_norm(trait_major_nlog.to_numpy()),
        xrot=0, show_stars=True, ylabels=True,
    )
    cb = fig.colorbar(im, ax=ax, fraction=0.038, pad=0.02)
    cb.set_label("-log10(p_cauchy)", rotation=270, labelpad=18)
    cb.outline.set_visible(False)
    fig.tight_layout()
    fig.savefig(major_dir / "trait_major_nlog10_p_cauchy_heatmap.png", dpi=FIG_DPI, bbox_inches="tight")
    fig.savefig(major_dir / "trait_major_nlog10_p_cauchy_heatmap.pdf", dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)

    # 5. Trait spatial maps in the HPF region of slice 2.
    hpf_dir = OUTPUT_DIR / "05_slice2_HPF_spatial_maps"
    adata2_hpf = adata2[adata2.obs["region"].astype(str) == "HPF"].copy()
    for trait in TRAIT_KEYS:
        tag = sanitize_key(trait)
        plot_spatial(
            adata2_hpf,
            trait,
            hpf_dir / f"slice2_HPF.{tag}.png",
            hpf_dir / f"slice2_HPF.{tag}.pdf",
            cmap_bwr,
        )
    plot_spatial_overview(
        adata2_hpf,
        TRAIT_KEYS,
        hpf_dir / "slice2_HPF_overview.png",
        hpf_dir / "slice2_HPF_overview.pdf",
        cmap_bwr,
    )

    # 6-7. Lee's L and bivariate Moran's I compared with cell-type enrichment.
    lee_mat, moran_bv_mat = build_spatial_trait_mats(
        adata1, adata2, TRAIT_KEYS_NO_TOTAL, k=K_SPATIAL
    )
    lee_dir = OUTPUT_DIR / "06_lee_with_major_heatmap"
    moran_dir = OUTPUT_DIR / "07_moran_bv_with_major_heatmap"
    save_df(lee_mat, lee_dir / "lee_spatial_association.tsv")
    save_df(moran_bv_mat, moran_dir / "moran_bv_spatial_association.tsv")

    plot_dual_heatmap(
        left_df=lee_mat,
        right_df=trait_major_nlog,
        right_p_df=trait_major_p,
        left_label="Lee's L",
        right_label="-log10(p_cauchy)",
        left_cmap=cmap_pearson,
        right_cmap=cmap_pearson,
        out_png=lee_dir / "lee_with_major_heatmap.png",
        out_pdf=lee_dir / "lee_with_major_heatmap.pdf",
    )
    plot_dual_heatmap(
        left_df=moran_bv_mat,
        right_df=trait_major_nlog,
        right_p_df=trait_major_p,
        left_label="Bivariate Moran's I",
        right_label="-log10(p_cauchy)",
        left_cmap=cmap_pearson,
        right_cmap=cmap_pearson,
        out_png=moran_dir / "moran_bv_with_major_heatmap.png",
        out_pdf=moran_dir / "moran_bv_with_major_heatmap.pdf",
    )

    # 8. Anatomical region maps and shared legend.
    region_map_dir = OUTPUT_DIR / "08_region_maps"
    region_palette = build_region_palette(adata1, adata2)
    plot_region_pair(
        adata1,
        adata2,
        region_palette,
        region_map_dir / "two_slices_region_map.png",
        region_map_dir / "two_slices_region_map.pdf",
    )
    plot_region_legend(
        region_palette,
        adata1,
        adata2,
        region_map_dir / "region_legend_no_root.png",
        region_map_dir / "region_legend_no_root.pdf",
    )

    # 9. Region and major-cell-type enrichment heatmaps in a shared layout.
    region_p_no_total = region_p[
        [TRAIT_LABELS[trait] for trait in TRAIT_KEYS_NO_TOTAL]
    ].copy()
    region_nlog_no_total = region_nlog[
        [TRAIT_LABELS[trait] for trait in TRAIT_KEYS_NO_TOTAL]
    ].copy()
    region_p_t = region_p_no_total.T.reindex(index=trait_major_p.index)
    region_nlog_t = region_nlog_no_total.T.reindex(index=trait_major_nlog.index)

    transposed_dir = OUTPUT_DIR / "09_region_transposed_with_major_heatmap"
    save_df(region_p_t, transposed_dir / "region_trait_transposed_p_cauchy.tsv")
    save_df(region_nlog_t, transposed_dir / "region_trait_transposed_nlog10_p_cauchy.tsv")
    plot_dual_heatmap_transposed(
        left_df=region_nlog_t,
        left_p_df=region_p_t,
        right_df=trait_major_nlog,
        right_p_df=trait_major_p,
        left_label="-log10(p_cauchy)",
        right_label="-log10(p_cauchy)",
        left_cmap=cmap_bwr,
        right_cmap=cmap_pearson,
        out_png=transposed_dir / "region_transposed_with_major_heatmap.png",
        out_pdf=transposed_dir / "region_transposed_with_major_heatmap.pdf",
    )


if __name__ == "__main__":
    run_analysis()