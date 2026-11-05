"""Layer 1b: where each enzyme is expressed in skin and how disease/cytokines move it.

Contrasts (all public):
  AD_LS_vs_HC, AD_NL_vs_HC, AD_LS_vs_NL, PSO_LS_vs_HC  GSE121212 RNA-seq (Tsoi 2019)
  KC_IL4, KC_IL13, KC_IL17A, KC_IFNa, KC_IFNg, KC_TNF GSE36287 primary keratinocytes
  ACNE_LS_vs_NL                                       GSE53795
Cell-type localisation: Human Protein Atlas skin scRNA-seq clusters.

Statistics are deliberately simple and transparent: log2 CPM (RNA-seq) or RMA
log2 (arrays), Welch or paired t-tests, Benjamini-Hochberg FDR.
"""
from __future__ import annotations

import gzip

import numpy as np
import pandas as pd
from scipy import stats

from .data import RAW, read_series_matrix


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    n = len(pv)
    if n == 0:
        return out
    order = np.argsort(pv)
    ranked = pv[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    res = np.empty(n)
    res[order] = np.clip(ranked, 0, 1)
    out[ok] = res
    return out


def contrast(x: pd.DataFrame, a: list[str], b: list[str], paired: bool = False, name: str = "") -> pd.DataFrame:
    """log2FC (a - b), t-test p and BH q for every row of x (genes x samples, log scale)."""
    missing = [c for c in (*a, *b) if c not in x.columns]
    if missing:
        raise ValueError(f"contrast {name or '(unnamed)'}: samples not in data: {missing}")
    if len(a) < 2 or len(b) < 2:
        raise ValueError(f"contrast {name or '(unnamed)'}: need >=2 samples per group for a t-test; "
                         f"got {len(a)} vs {len(b)}")
    if paired and len(a) != len(b):
        raise ValueError(f"contrast {name or '(unnamed)'}: paired test needs equal group sizes; "
                         f"got {len(a)} vs {len(b)}")
    xa, xb = x[a].to_numpy(float), x[b].to_numpy(float)
    if paired:
        t, p = stats.ttest_rel(xa, xb, axis=1)
    else:
        t, p = stats.ttest_ind(xa, xb, axis=1, equal_var=False)
    lfc = np.nanmean(xa, 1) - np.nanmean(xb, 1)
    df = pd.DataFrame({"lfc": lfc, "p": p, "q": bh(p)}, index=x.index)
    return df.add_prefix(f"{name}_") if name else df


# ---------------------------------------------------------------- GSE121212
def ad_psoriasis() -> pd.DataFrame:
    counts = pd.read_csv(RAW / "geo" / "GSE121212_readcount.txt.gz", sep="\t", index_col=0)
    counts = counts.groupby(level=0).sum()
    lib = counts.sum(axis=0)
    logcpm = np.log2(counts / lib * 1e6 + 1)
    logcpm = logcpm[(counts / lib * 1e6 >= 1).sum(axis=1) >= 10]  # expressed in >=10 samples

    cols = pd.Series(logcpm.columns)
    ad_ls = cols[cols.str.match(r"AD_.*_lesional")].tolist()
    ad_nl = cols[cols.str.match(r"AD_.*_non-lesional")].tolist()
    pso_ls = cols[cols.str.match(r"PSO_.*_lesional")].tolist()
    hc = cols[cols.str.match(r"CTRL_")].tolist()
    # Paired AD LS vs NL on patients with both biopsies
    pid = lambda s: s.rsplit("_", 1)[0]
    both = sorted(set(map(pid, ad_ls)) & set(map(pid, ad_nl)))
    ls_p = [f"{p}_lesional" for p in both]
    nl_p = [f"{p}_non-lesional" for p in both]

    out = pd.concat([
        contrast(logcpm, ad_ls, hc, name="AD_LS_vs_HC"),
        contrast(logcpm, ad_nl, hc, name="AD_NL_vs_HC"),
        contrast(logcpm, ls_p, nl_p, paired=True, name="AD_LS_vs_NL"),
        contrast(logcpm, pso_ls, hc, name="PSO_LS_vs_HC"),
    ], axis=1)
    out.attrs["n"] = dict(AD_LS=len(ad_ls), AD_NL=len(ad_nl), PSO_LS=len(pso_ls), HC=len(hc), AD_pairs=len(both))
    out["skin_log2cpm_HC"] = logcpm[hc].mean(axis=1)
    return out


# ---------------------------------------------------------------- GPL570 arrays
def _gpl570_gene_map() -> pd.Series:
    with gzip.open(RAW / "geo" / "GPL570.annot.gz", "rt", errors="replace") as fh:
        lines = [l for l in fh if not l.startswith(("!", "#", "^"))]
    from io import StringIO

    a = pd.read_csv(StringIO("".join(lines)), sep="\t", usecols=["ID", "Gene symbol"], dtype=str)
    a = a.dropna()
    a = a[~a["Gene symbol"].str.contains("///")]
    return a.set_index("ID")["Gene symbol"]


def _collapse_probes(expr: pd.DataFrame, gmap: pd.Series) -> pd.DataFrame:
    e = expr.loc[expr.index.intersection(gmap.index)]
    if e.max().max() > 50:  # not yet log-transformed
        e = np.log2(e.clip(lower=1))
    e = e.assign(_g=gmap.reindex(e.index).values, _m=e.mean(axis=1))
    # keep the most highly expressed probe per gene
    e = e.sort_values("_m", ascending=False).drop_duplicates("_g").set_index("_g").drop(columns="_m")
    e.index.name = "gene"
    return e


def keratinocyte_cytokines() -> pd.DataFrame:
    meta, expr = read_series_matrix(RAW / "geo" / "GSE36287_series_matrix.txt.gz")
    x = _collapse_probes(expr, _gpl570_gene_map())
    title = dict(zip(meta["geo_accession"], meta["title"]))
    x.columns = [title[c].replace("Keratinocytes_", "") for c in x.columns]
    outs = []
    for cyt in ["IL4", "IL13", "IFNa", "IFNg", "TNF"]:
        trt = [f"{cyt}_subject{i}" for i in (1, 2, 3)]
        ctl = [f"untreated_subject{i}" for i in (1, 2, 3)]
        outs.append(contrast(x, trt, ctl, paired=True, name=f"KC_{cyt}"))
    # IL-17A was run as a separate batch with its own untreated controls
    outs.append(contrast(x, [f"IL17A_replicate{i}" for i in (1, 2, 3)],
                         [f"untreated_replicate{i}" for i in (1, 2, 3)], paired=True, name="KC_IL17A"))
    out = pd.concat(outs, axis=1)
    out["KC_log2_untreated"] = x[[f"untreated_subject{i}" for i in (1, 2, 3)]].mean(axis=1)
    return out


def acne() -> pd.DataFrame:
    meta, expr = read_series_matrix(RAW / "geo" / "GSE53795_series_matrix.txt.gz")
    x = _collapse_probes(expr, _gpl570_gene_map())
    title = dict(zip(meta["geo_accession"], meta["title"]))
    x.columns = [title[c] for c in x.columns]
    pts = sorted({c.split("_Pt")[1] for c in x.columns})
    ls = [f"Acne lesion_Pt{p}" for p in pts if f"Acne lesion_Pt{p}" in x and f"Acne non lesional_Pt{p}" in x]
    nl = [c.replace("Acne lesion", "Acne non lesional") for c in ls]
    return contrast(x, ls, nl, paired=True, name="ACNE_LS_vs_NL")


# ---------------------------------------------------------------- HPA skin scRNA
KC_TYPES = ["basal keratinocytes", "suprabasal keratinocytes"]


def skin_cell_types() -> pd.DataFrame:
    h = pd.read_csv(RAW / "hpa" / "skin_clusters_ntpm.tsv", sep="\t")
    # cell-count-weighted mean nCPM per cell type
    cl = pd.read_csv(RAW / "hpa" / "rna_single_cell_clusters.tsv", sep="\t")
    cl = cl[cl["Tissue"] == "skin"][["Cluster", "Cell count"]]
    h = h.merge(cl, on="Cluster")
    h["w"] = h["nCPM"] * h["Cell count"]
    g = h.groupby(["Gene name", "Cell type"])[["w", "Cell count"]].sum()
    ct = (g["w"] / g["Cell count"]).unstack().fillna(0)
    ct.index.name = "gene"
    out = pd.DataFrame(index=ct.index)
    out["top_skin_cell_type"] = ct.idxmax(axis=1)
    out["top_skin_nCPM"] = ct.max(axis=1)
    for t in KC_TYPES + ["melanocytes", "fibroblasts", "sebaceous gland cells"]:
        if t in ct:
            out[f"nCPM_{t.replace(' ', '_')}"] = ct[t]
    kc = ct[[t for t in KC_TYPES if t in ct]].max(axis=1)
    out["KC_nCPM"] = kc
    # keratinocyte specificity: share of the gene's expression that sits in KCs
    out["KC_specificity"] = kc / ct.sum(axis=1).replace(0, np.nan)
    # differentiation direction: suprabasal vs basal
    if set(KC_TYPES) <= set(ct):
        out["suprabasal_vs_basal_log2"] = np.log2((ct["suprabasal keratinocytes"] + 1) / (ct["basal keratinocytes"] + 1))
    return out


def build_expression() -> pd.DataFrame:
    parts = [ad_psoriasis(), keratinocyte_cytokines(), acne(), skin_cell_types()]
    return pd.concat(parts, axis=1)
