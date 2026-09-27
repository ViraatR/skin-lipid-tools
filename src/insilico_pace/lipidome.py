"""Layer 2: infer which enzyme activities change from the lipidome.

Primary method: product/substrate ratios ("activity proxies"). Each proxy is a
log2 ratio of summed species on either side of a known enzymatic step, computed
per animal, then tested for a genotype effect with age as a covariate. This is
the same logic LINEX2 uses on individual reaction edges, collapsed to the
class/chain level so it stays interpretable with n = 3-4 animals per group.
`linex_enrichment()` runs LINEX2 itself for comparison.

Dataset: Metabolomics Workbench ST002195, mouse Tmem79 (Spade) dermatitis
model vs WT dorsal skin at P0, 4, 8 and 10 weeks (RPLC-MS, MS-DIAL names).
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Callable

import numpy as np
import pandas as pd
from scipy import stats

from lipidname import Lipid, parse

from .data import load_mwb
from .transcriptome import bh

SampleFilter = Callable[[Lipid], bool]


@dataclass
class Proxy:
    name: str
    numerator: tuple[str | tuple[str, ...], SampleFilter | None]
    denominator: tuple[str | tuple[str, ...], SampleFilter | None]
    up: tuple[str, ...]      # enzymes whose activity raises the ratio
    down: tuple[str, ...] = ()  # enzymes whose activity lowers it
    note: str = ""


def _acyl(lo: int, hi: int) -> SampleFilter:
    return lambda l: l.acyl_c is not None and lo <= l.acyl_c <= hi


def _chains(pairs: set[tuple[int, int]]) -> SampleFilter:
    """Weight = number of acyl chains in `pairs` (needs molecular-species annotation)."""
    return lambda l: sum(ch in pairs for ch in l.chains) if len(l.chains) > 1 else 0


GLYCEROLIPIDS = ("PC", "PE", "PS", "PI", "PG", "TG", "DG")


def _fa_len(lo: int, hi: int) -> SampleFilter:
    return lambda l: l.total_c is not None and lo <= l.total_c <= hi


PLA2 = ("PLA2G4A", "PLA2G4B", "PLA2G4D", "PLA2G4E", "PLA2G4F", "PLA2G6", "PNPLA8", "PLA2G2F", "PLA2G10", "PLA2G15", "LYPLA1")
LPCAT = ("LPCAT1", "LPCAT2", "LPCAT3", "LPCAT4", "MBOAT7", "MBOAT2")

PROXIES = [
    Proxy("SM_hydrolysis", ("Cer_NS", None), ("SM", None),
          up=("SMPD1", "SMPD2", "SMPD3", "SMPD4", "ENPP7"), down=("SGMS1", "SGMS2"),
          note="Cer[NS]/SM: sphingomyelinase vs SM synthase"),
    Proxy("Cer_glucosylation", ("HexCer_NS", None), ("Cer_NS", None),
          up=("UGCG",), down=("GBA1", "GBA2"),
          note="GlcCer/Cer: the SC ceramide pool is released from GlcCer by GBA1"),
    Proxy("DHCer_desaturation", ("Cer_NS", None), ("Cer_NDS", None), up=("DEGS1",),
          note="Cer[NS]/Cer[NDS]: sphinganine -> sphingosine backbone"),
    Proxy("Phyto_hydroxylation", ("Cer_NP", None), ("Cer_NDS", None), up=("DEGS2",),
          note="Cer[NP]/Cer[NDS]: C4-hydroxylation to phytoceramide"),
    Proxy("Alpha_hydroxylation", (("Cer_AS", "Cer_ADS", "Cer_AP"), None), (("Cer_NS", "Cer_NDS", "Cer_NP"), None),
          up=("FA2H",), note="alpha-OH / non-OH ceramides"),
    Proxy("Acylceramide", ("Cer_EOS", None), ("Cer_NS", _acyl(28, 40)),
          up=("PNPLA1", "CYP4F22", "ALOX12B", "ALOXE3", "SDR9C7", "ABHD5"),
          note="Cer[EOS]/ultra-long Cer[NS]: omega-O-acylceramide pathway"),
    Proxy("ULC_ceramide", ("Cer_NS", _acyl(26, 40)), ("Cer_NS", _acyl(16, 24)),
          up=("CERS3", "ELOVL4"), note="acyl >=C26 vs <=C24 in Cer[NS]"),
    Proxy("VLC_ceramide", ("Cer_NS", _acyl(22, 24)), ("Cer_NS", _acyl(16, 18)),
          up=("CERS2", "ELOVL1"), down=("CERS5", "CERS6"), note="C22-24 vs C16-18 Cer[NS]"),
    Proxy("C18_C16_chains", (GLYCEROLIPIDS, _chains({(18, 0), (18, 1)})), (GLYCEROLIPIDS, _chains({(16, 0), (16, 1)})),
          up=("ELOVL6",), note="18:0+18:1 vs 16:0+16:1 acyl chains in glycerolipids (free C16/C18 FA not reported)"),
    Proxy("FA_VLC", ("FA", _fa_len(20, 36)), ("FA", _fa_len(14, 18)),
          up=("ELOVL1", "ELOVL3", "ELOVL7"), note=">=C20 vs C14-18 free fatty acids"),
    Proxy("SCD_index", (GLYCEROLIPIDS, _chains({(16, 1), (18, 1)})), (GLYCEROLIPIDS, _chains({(16, 0), (18, 0)})),
          up=("SCD", "SCD5"), note="(16:1+18:1)/(16:0+18:0) acyl chains in glycerolipids"),
    Proxy("PC_deacylation", ("LPC", None), ("PC", None), up=PLA2, down=LPCAT, note="LPC/PC"),
    Proxy("PE_deacylation", ("LPE", None), ("PE", None), up=PLA2, down=LPCAT, note="LPE/PE"),
    Proxy("PE_methylation", ("PC", None), ("PE", None), up=("PEMT",), note="PC/PE"),
    Proxy("PS_decarboxylation", ("PE", None), ("PS", None), up=("PISD",), down=("PTDSS2",), note="PE/PS"),
    Proxy("Plasmalogen", ("PE-P", None), ("PE", None), up=("FAR1", "GNPAT", "AGPS", "PEDS1"), note="PE-P/PE"),
    Proxy("TG_lipolysis", ("DG", None), ("TG", None), up=("PNPLA2", "LIPE"), down=("DGAT1", "DGAT2"),
          note="DG/TG"),
]


def load_spade() -> tuple[pd.DataFrame, pd.DataFrame]:
    """-> (features with parsed Lipid, animal x feature linear intensities)."""
    feats, mat = load_mwb("ST002195")
    feats["lipid"] = feats["metabolite_name"].map(parse)
    feats = feats[feats["lipid"].notna()].copy()
    feats["mode"] = feats["analysis_summary"].str.contains("NEG").map({True: "neg", False: "pos"})
    feats["cls"] = feats["lipid"].map(lambda l: l.lipid_class)

    # Each feature is only measured in its own polarity's injections; stitch per animal
    per_mode = {}
    for mode in ("neg", "pos"):
        ids = feats.index[feats["mode"] == mode]
        m = mat.loc[mat.index.str.endswith(f"_{mode}"), ids]
        m.index = m.index.str.replace(f"_{mode}$", "", regex=True)
        per_mode[mode] = m
    x = pd.concat(per_mode.values(), axis=1)
    # Half-minimum imputation per feature for non-detects (standard for sparse lipidomics)
    x = x.apply(lambda c: c.fillna(c[c > 0].min() / 2 if (c > 0).any() else 0)).clip(lower=0)
    return feats, x


def sample_design(index: pd.Index) -> pd.DataFrame:
    parts = index.str.extract(r"^(?P<id>\d+)_(?P<genotype>sp|WT)_(?P<age>\w+)$")
    parts.index = index
    parts["spade"] = (parts["genotype"] == "sp").astype(int)
    return parts


def _class_sum(feats, x, classes, filt) -> pd.Series:
    classes = (classes,) if isinstance(classes, str) else classes
    sel = feats[feats["cls"].isin(classes)]
    # a filter returns bool (keep species) or a number (weight, e.g. matching chain count)
    w = sel["lipid"].map(filt).astype(float) if filt is not None else pd.Series(1.0, index=sel.index)
    sel = sel[w > 0]
    if sel.empty:
        return pd.Series(np.nan, index=x.index)
    # use one polarity per class (the one with more species) to avoid double counting
    best = sel.groupby(["cls", "mode"]).size().reset_index().sort_values(0).drop_duplicates("cls", keep="last")
    chosen = set(zip(best["cls"], best["mode"]))
    keep = sel.index[[(c, m) in chosen for c, m in zip(sel["cls"], sel["mode"])]]
    return (x[keep] * w[keep]).sum(axis=1)


def compute_proxies(feats: pd.DataFrame, x: pd.DataFrame) -> pd.DataFrame:
    out = {}
    for p in PROXIES:
        num = _class_sum(feats, x, *p.numerator)
        den = _class_sum(feats, x, *p.denominator)
        out[p.name] = np.log2(num / den)
    return pd.DataFrame(out)


def _ols_genotype(y: pd.Series, d: pd.DataFrame) -> tuple[float, float, float]:
    """y ~ const + spade + age dummies; returns (spade coef, t, p). NaNs if not estimable."""
    X = pd.get_dummies(d["age"], drop_first=True, dtype=float)
    X.insert(0, "spade", d["spade"].astype(float))
    X.insert(0, "const", 1.0)
    X, yy = X.to_numpy(), y.to_numpy(float)
    dof = len(yy) - X.shape[1]
    if dof < 2 or np.linalg.matrix_rank(X) < X.shape[1]:
        return np.nan, np.nan, np.nan
    beta, *_ = np.linalg.lstsq(X, yy, rcond=None)
    resid = yy - X @ beta
    se = np.sqrt(resid @ resid / dof * np.linalg.inv(X.T @ X)[1, 1])
    if se == 0:
        return beta[1], 0.0, 1.0
    t = beta[1] / se
    return beta[1], t, 2 * stats.t.sf(abs(t), dof)


def genotype_effects(proxies: pd.DataFrame, design: pd.DataFrame, ages=("4w", "8w", "10w")) -> pd.DataFrame:
    """OLS per proxy: log2 ratio ~ spade + age; reports the spade coefficient.

    P0 is excluded by default: dermatitis in Spade mice develops after weaning,
    so P0 serves as the pre-disease reference and is reported separately.
    """
    d = design[design["age"].isin(ages)]
    rows = []
    for name in proxies:
        y = proxies.loc[d.index, name]
        ok = y.notna() & np.isfinite(y)
        beta, t, p = _ols_genotype(y[ok], d[ok])
        per_age = {f"lfc_{a}": proxies.loc[design.index[(design.age == a) & (design.spade == 1)], name].mean()
                   - proxies.loc[design.index[(design.age == a) & (design.spade == 0)], name].mean()
                   for a in ["P0", *ages]}
        rows.append({"proxy": name, "lfc_spade": beta, "t": t, "p": p, **per_age})
    res = pd.DataFrame(rows).set_index("proxy")
    res["q"] = bh(res["p"].to_numpy())
    meta = pd.DataFrame([{"proxy": p.name, "note": p.note, "up": ",".join(p.up), "down": ",".join(p.down)} for p in PROXIES]).set_index("proxy")
    return res.join(meta)


def enzyme_lipid_evidence(effects: pd.DataFrame) -> pd.DataFrame:
    """Project proxy effects onto enzymes: signed t (+ = inferred activity up in dermatitis)."""
    rows = []
    for proxy, r in effects.iterrows():
        for sign, genes in ((1, r["up"]), (-1, r["down"])):
            for g in filter(None, genes.split(",")):
                rows.append({"gene": g, "proxy": proxy, "activity_t": sign * r["t"], "q": r["q"],
                             "activity_lfc": sign * r["lfc_spade"]})
    ev = pd.DataFrame(rows)
    # keep the strongest proxy per enzyme
    ev = ev.loc[ev.groupby("gene")["activity_t"].apply(lambda s: s.abs().idxmax())]
    return ev.set_index("gene").rename(columns={"proxy": "lipid_proxy", "q": "lipid_q"})


def species_effects(feats: pd.DataFrame, x: pd.DataFrame, design: pd.DataFrame, ages=("4w", "8w", "10w")) -> pd.DataFrame:
    """Per-species Spade vs WT (age-adjusted) for the supplementary table / LINEX2."""
    lx = np.log2(x.replace(0, np.nan)).fillna(0)
    d = design[design["age"].isin(ages)]
    rows = []
    for fid in feats.index:
        beta, t, p = _ols_genotype(lx.loc[d.index, fid], d)
        rows.append({"feature": fid, "name": feats.at[fid, "metabolite_name"], "species": feats.at[fid, "lipid"].key,
                     "class": feats.at[fid, "cls"], "lfc_spade": beta, "t": t, "p": p})
    out = pd.DataFrame(rows).set_index("feature")
    out["q"] = bh(out["p"].to_numpy())
    return out.sort_values("p")


# -------------------------------------------------------------------- LINEX2
def _to_linex(l: Lipid) -> str | None:
    """Sum-composition name in the 'PC(34:1)' / 'Cer(d42:1)' form LINEX2 parses
    without LipidLynxX. Hydroxy-acyl (A/B) and acyl- (EO) ceramides have no
    distinct LINEX2 class, so they are left out rather than mislabelled."""
    cls, c, db = l.lipid_class, l.total_c, l.total_db
    if cls.startswith("Cer_"):
        sub = cls[4:]
        if sub.startswith(("A", "B", "EO", "O")):
            return None
        return f"Cer({'t' if sub.endswith('P') else 'd'}{c}:{db})"
    if cls.startswith("HexCer_"):
        return f"HexCer(d{c}:{db})"
    if cls == "SM":
        return f"SM(d{c}:{db})"
    if cls in ("PE-P", "PC-P"):  # plasmalogen P-x:y == O-x:(y+1)
        return f"{cls[:2]}(O-{c}:{db + 1})"
    if cls in ("PC", "PE", "PS", "PI", "PG", "LPC", "LPE", "LPS", "LPI", "LPG", "TG", "DG", "CE", "PA", "LPA"):
        return f"{cls}({c}:{db})"
    return None


def linex_enrichment(feats: pd.DataFrame, x: pd.DataFrame, design: pd.DataFrame, ages=("8w", "10w"), **kw):
    """Run LINEX2 reaction-network enrichment (Spade vs WT). Returns (enzymes, score, p)."""
    import contextlib
    import io

    from linex2 import LipidNetwork

    d = design[design["age"].isin(ages)]
    names = feats["lipid"].map(_to_linex)
    keep = names.dropna()
    data = np.log2(x.loc[d.index, keep.index].replace(0, np.nan)).fillna(0)
    data.columns = keep.values
    data = data.T.groupby(level=0).max().T  # collapse duplicate species across modes
    groups = d["genotype"].map({"sp": "Spade", "WT": "WT"})
    with contextlib.redirect_stdout(io.StringIO()):
        net = LipidNetwork(data, groups=groups, lipids_are_columns=True, organism="MMU",
                           database=("Rhea", "Reactome"), allow_molspec_fails=True)
        net.compute_native_network()
        sub, scores, pvals, _ = net.analysis_reaction_network_enrichment(
            "WT", "Spade", plot_scores=False, repetitions=kw.pop("repetitions", 3),
            p_value_iterations=kw.pop("p_value_iterations", 2000), **kw)
    return net, sub, scores, pvals


def linex_enzymes(subnetworks: dict, pvals: dict, rhea_to_genes: dict[str, set[str]] | None = None) -> pd.DataFrame:
    """Gene -> number of enriched-subnetwork reaction edges it annotates (+ component p).

    LINEX2 edges carry gene names only sometimes; Rhea IDs are mapped to human
    genes through the UniProt catalog (`rhea_to_genes`) to fill the gaps.
    """
    import re as _re

    rhea_to_genes = rhea_to_genes or {}
    rows = []
    for comp, g in subnetworks.items():
        for u, v, e in g.edges(data=True):
            genes = {x.strip().upper() for x in str(e.get("enzyme_gene_name", "")).split(",") if x.strip()}
            for rid in _re.findall(r"RHEA:(\d+)", str(e.get("enzyme_id", ""))):
                genes |= rhea_to_genes.get(rid, set())
            for gene in genes:
                rows.append({"gene": gene, "component": comp, "reaction": e.get("reaction_type"),
                             "edge": f"{u} -> {v}", "linex_p": pvals.get(comp)})
    if not rows:
        return pd.DataFrame(columns=["linex_edges", "linex_p", "linex_reactions"])
    df = pd.DataFrame(rows)
    return df.groupby("gene").agg(linex_edges=("edge", "nunique"), linex_p=("linex_p", "min"),
                                  linex_reactions=("reaction", lambda s: ",".join(sorted(set(s)))))
