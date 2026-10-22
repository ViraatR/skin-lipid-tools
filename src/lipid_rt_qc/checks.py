"""Annotation and data-quality checks. Each returns a per-feature DataFrame.

1. rt_outliers        observed RT disagrees with the dataset's own RT model
2. elution_order      within a class, RT must rise with carbons and fall with
                      double bonds (reverse phase). Includes the Cer vs DHCer rule
                      from the lab's Skyline notes: Cer[NS] x:y elutes before
                      Cer[NDS] x:(y-1).
3. isobars            other plausible annotations at the same m/z, and whether
                      the RT model prefers one of them
4. duplicates         one annotation at several RTs: split peak, positional
                      isomer, or a monomethyl branched-chain (mmBCFA) isomer,
                      which elutes slightly earlier than its straight-chain form
5. blank / QC CV      the filtering usually done by hand
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from lipidname import parse

from .io import Dataset
from .masses import ADDUCTS, NEGATIVE, POSITIVE, enumerate_isobars, mz as calc_mz
from .model import RTModel, expected_chains, lipid_features, model_class, robust_z


def annotate(ds: Dataset) -> pd.DataFrame:
    f = ds.features.copy()
    f["lipid"] = f["name"].map(parse)
    f["lipid_class"] = f["lipid"].map(lambda l: model_class(l) if l else None)
    f["species"] = f["lipid"].map(lambda l: l.key if l else None)
    # MS-DIAL puts the adduct in the name when there is no separate column
    missing = f["adduct"].fillna("").eq("")
    f.loc[missing, "adduct"] = f.loc[missing, "lipid"].map(lambda l: (l.adduct if l else "") or "")
    return f


# ------------------------------------------------------------------ 1. RT model
def rt_outliers(f: pd.DataFrame, z_thresh: float = 3.0, min_dt: float = 0.2, model: RTModel | None = None):
    X = lipid_features(f["lipid"].dropna())
    X = X[f.loc[X.index, "rt"].notna()]
    if len(X) < 2:
        raise ValueError(f"RT model needs at least 2 features with a parseable lipid name and an RT; got {len(X)}")
    y = f.loc[X.index, "rt"]
    model = model or RTModel()
    pred = pd.Series(model.oof_predict(X, y), index=X.index)
    out = pd.DataFrame(index=f.index)
    out["rt_pred"] = pred
    out["rt_resid"] = f["rt"] - pred
    out["rt_z"] = robust_z(out["rt_resid"], f["lipid_class"].fillna("?"))
    out["flag_rt"] = (out["rt_z"].abs() > z_thresh) & (out["rt_resid"].abs() > min_dt)
    return out, model, X


# ------------------------------------------------------------------ 2. elution order
def elution_order(f: pd.DataFrame, rt: pd.DataFrame, tol: float = 0.1) -> pd.DataFrame:
    """Reverse-phase ordering rules, evaluated on one representative RT per name.

    Every violating pair is blamed on the member that fits the RT model worse, so a
    single mis-annotation does not also flag the correct species it collides with.
    """
    ok = f.dropna(subset=["lipid", "rt"]).join(rt["rt_z"])
    ok["az"] = ok["rt_z"].abs().fillna(0)
    blame = pd.Series(0, index=f.index)
    notes: dict[str, list[str]] = {i: [] for i in f.index}

    def charge(a, b, why_a, why_b):
        if a.az >= b.az:
            blame[a.Index] += 1; notes[a.Index].append(why_a)
        else:
            blame[b.Index] += 1; notes[b.Index].append(why_b)

    ok = ok.assign(cls=ok["lipid"].map(model_class), c=ok["lipid"].map(lambda l: l.total_c),
                   d=ok["lipid"].map(lambda l: l.total_db))
    for cls, g in ok.groupby("cls"):
        rows = list(g.itertuples())
        for a in rows:
            for b in rows:
                if a.Index >= b.Index:
                    continue
                lo, hi = (a, b) if a.c < b.c else (b, a)
                if a.d == b.d and a.c != b.c and hi.rt < lo.rt - tol:
                    charge(hi, lo, f"elutes before {cls} {lo.c}:{lo.d} despite +{hi.c - lo.c}C",
                           f"elutes after {cls} {hi.c}:{hi.d} despite -{hi.c - lo.c}C")
                lo, hi = (a, b) if a.d < b.d else (b, a)
                if a.c == b.c and a.d != b.d and hi.rt > lo.rt + tol:
                    charge(hi, lo, f"elutes after {cls} {lo.c}:{lo.d} despite +{hi.d - lo.d}DB",
                           f"elutes before {cls} {hi.c}:{hi.d} despite -{hi.d - lo.d}DB")
    # Cer vs DHCer (lab Skyline note): Cer[NS] C:D must elute before Cer[NDS] C:(D-1)
    ns = ok[ok["cls"] == "Cer_NS"]
    nds = ok[ok["cls"] == "Cer_NDS"]
    for a in ns.itertuples():
        for b in nds[(nds["c"] == a.c) & (nds["d"] == a.d - 1)].itertuples():
            if a.rt > b.rt + tol:
                charge(a, b, f"Cer elutes after its DHCer (Cer_NDS {b.c}:{b.d})",
                       f"DHCer elutes before its Cer (Cer_NS {a.c}:{a.d})")
    out = pd.DataFrame({"order_violations": blame,
                        "order_notes": pd.Series({k: "; ".join(dict.fromkeys(v))[:300] for k, v in notes.items()})})
    out["flag_order"] = out["order_violations"] >= 2
    return out


# ------------------------------------------------------------------ 3. isobars
def _polarity(adduct: str) -> set[str]:
    a = (adduct or "").strip()
    if a.endswith("-"):
        return NEGATIVE
    if a.endswith("+"):
        return POSITIVE
    return POSITIVE | NEGATIVE


def isobars(f: pd.DataFrame, rt: pd.DataFrame, model: RTModel, X: pd.DataFrame, ppm: float = 5.0,
            classes_seen_only: bool = True, z_thresh: float = 3.0) -> pd.DataFrame:
    """Alternative annotations within `ppm`; score each alternative with the RT model."""
    from .masses import CLASS_RULES
    classes = sorted(c for c in X["lipid_class"].unique() if c in CLASS_RULES) if classes_seen_only else None
    out = pd.DataFrame(index=f.index, columns=["n_isobars", "isobars", "best_alt", "best_alt_z", "flag_isobar_rt"])
    resid_sd = 1.4826 * np.nanmedian(np.abs(rt["rt_resid"] - np.nanmedian(rt["rt_resid"])))
    for fid, r in f.iterrows():
        l = r["lipid"]
        if l is None or pd.isna(r["rt"]):
            continue
        target = r["mz"] if pd.notna(r["mz"]) else (calc_mz(l.lipid_class, l.total_c, l.total_db, r["adduct"])
                                                     if r["adduct"] in ADDUCTS else None)
        if target is None:
            continue
        alts = [a for a in enumerate_isobars(target, _polarity(r["adduct"]), ppm, classes)
                if not (a[0] == l.lipid_class and a[1] == l.total_c and a[2] == l.total_db)]
        # same-formula isomers across Cer subclasses (NS/NDS) are tagged but are not "errors"
        out.at[fid, "n_isobars"] = len(alts)
        if not alts:
            continue
        cand = _alt_features(pd.DataFrame(alts, columns=["lipid_class", "total_c", "total_db", "adduct", "ppm"]), X)
        cand["rt_pred"] = model.predict(cand)
        cand["z"] = (r["rt"] - cand["rt_pred"]) / (resid_sd or 1)
        cand = cand.sort_values(by="z", key=abs)
        out.at[fid, "isobars"] = "; ".join(
            f"{c.lipid_class} {c.total_c}:{c.total_db} {c.adduct} ({c.ppm:+.1f} ppm, z={c.z:+.1f})" for c in cand.head(4).itertuples())
        b = cand.iloc[0]
        out.at[fid, "best_alt"] = f"{b.lipid_class} {b.total_c}:{b.total_db} {b.adduct}"
        out.at[fid, "best_alt_z"] = b.z
    cur_z = rt["rt_z"].abs()
    alt_z = pd.to_numeric(out["best_alt_z"], errors="coerce").abs()
    # flag when the current annotation fits RT poorly and an isobar fits well
    out["flag_isobar_rt"] = (cur_z > z_thresh) & (alt_z < 1.5)
    out["n_isobars"] = pd.to_numeric(out["n_isobars"], errors="coerce").fillna(0).astype(int)
    return out


def _alt_features(cand: pd.DataFrame, X: pd.DataFrame) -> pd.DataFrame:
    """Fill the RT-model features an alternative annotation would have.

    Classes seen in the run take the typical OH / ether / base-DB / chain values of
    their own annotated features, so e.g. a Cer_NP alternative gets the phyto base
    (no base double bond, extra OH) the model was trained on. Unseen classes fall
    back to rules derived from the class name.
    """
    cand = cand.copy()
    cols = ["total_oh", "ether", "base_db", "n_chains"]
    cand["total_oh"] = 0
    cand["ether"] = cand["lipid_class"].str.endswith(("-O", "-P")).astype(int)
    cand["base_db"] = np.where(cand["lipid_class"].str.startswith(("Cer", "SM", "HexCer")),
                               np.where(cand["lipid_class"].str.contains("DS"), 0, 1), -1)
    cand["n_chains"] = cand["lipid_class"].map(expected_chains)
    typical = X.groupby("lipid_class")[cols].agg(lambda v: v.mode().iloc[0])
    seen = cand["lipid_class"].isin(typical.index)
    if seen.any():
        cand.loc[seen, cols] = typical.loc[cand.loc[seen, "lipid_class"], cols].to_numpy()
    return cand


# ------------------------------------------------------------------ 4. duplicates / BCFA
def duplicates(f: pd.DataFrame, rt: pd.DataFrame, X: pd.DataFrame, model: RTModel,
               split_peak_dt: float = 0.1) -> pd.DataFrame:
    """Same species annotated at several RTs. Uses the class's per-carbon RT slope
    to express the gap in carbon units: mmBCFA (iso/anteiso) isomers elute about
    0.3-1 carbon-equivalents earlier than the straight-chain species."""
    out = pd.DataFrame(index=f.index, columns=["dup_group", "dup_n", "dup_dt", "dup_carbon_eq", "dup_call"])
    slope = _carbon_slope(X, f["rt"])
    key = f["name"].str.replace(r";\s*\[M.*$", "", regex=True).str.strip()
    for sp, g in f.assign(_key=key).dropna(subset=["lipid", "rt"]).groupby(["_key", "adduct"]):
        if len(g) < 2:
            continue
        g = g.sort_values("rt")
        latest = g["rt"].iloc[-1]
        cls = g["lipid_class"].iloc[0]
        s = slope.get(cls, np.nanmedian(list(slope.values())) if slope else np.nan)
        l0 = g["lipid"].iloc[0]
        # mmBCFAs are saturated; a species with no possible saturated chain can't carry one
        sat_possible = (not l0.chains) or any(d == 0 for _, d in l0.chains)
        for fid, rtv in g["rt"].items():
            dt = rtv - latest
            ceq = dt / s if s and not np.isnan(s) else np.nan
            if fid == g.index[-1]:
                call = "latest-eluting (straight-chain / main isomer)"
            elif abs(dt) <= split_peak_dt:
                call = "split peak? (same RT within tolerance)"
            elif l0.labelled:
                call = "labelled standard eluting twice: check integration"
            elif cls.endswith(("-O", "-P")):
                call = "ether isomer (plasmanyl O- vs plasmenyl P- share a formula)"
            elif -1.0 <= ceq <= -0.3 and sat_possible:
                call = "earlier by ~0.3-1 C-eq: branched-chain (mmBCFA) isomer candidate"
            else:
                call = "separate isomer (sn-position / double-bond position) or misannotation"
            out.loc[fid] = [f"{sp[0]} {sp[1]}".strip(), len(g), dt, ceq, call]
    out["flag_bcfa"] = out["dup_call"].fillna("").str.contains("branched")
    return out


def _carbon_slope(X: pd.DataFrame, rt: pd.Series) -> dict[str, float]:
    """Median RT increase per added carbon (at fixed DB) within each class."""
    s = {}
    d = X.join(rt.rename("rt")).dropna(subset=["rt"])
    for cls, g in d.groupby("lipid_class"):
        vals = []
        for db, gg in g.groupby("total_db"):
            if gg["total_c"].nunique() >= 2:
                p = np.polyfit(gg["total_c"], gg["rt"], 1)[0]
                if p > 0:
                    vals.append(p)
        if vals:
            s[cls] = float(np.median(vals))
    return s


# ------------------------------------------------------------------ 5. blank / QC
def blank_and_qc(ds: Dataset, fold: float = 3.0, max_cv: float = 30.0, min_detect: float | None = None) -> pd.DataFrame:
    out = pd.DataFrame(index=ds.features.index)
    if ds.intensity.empty or ds.samples.empty:
        return out
    x = ds.intensity.reindex(ds.features.index)
    typ = ds.samples["type"].reindex(x.columns)
    smp, blk, qc = (x.loc[:, typ == t] for t in ("sample", "blank", "qc"))
    if smp.shape[1]:
        out["sample_mean"] = smp.mean(axis=1)
        out["detect_rate"] = (smp.fillna(0) > 0).mean(axis=1)
        if min_detect is not None:
            out["flag_low_detect"] = out["detect_rate"] < min_detect
    if blk.shape[1] and smp.shape[1]:
        out["blank_mean"] = blk.mean(axis=1)
        # a feature never seen in samples has a NaN mean (e.g. Skyline replicates with no
        # peak); treat that as 0 so blank-only contaminants are still flagged
        out["sample_blank_ratio"] = out["sample_mean"].fillna(0) / out["blank_mean"].replace(0, np.nan)
        out["flag_blank"] = out["sample_blank_ratio"] < fold
    if qc.shape[1] >= 2:
        out["qc_cv"] = 100 * qc.std(axis=1, ddof=1) / qc.mean(axis=1)
        out["flag_qc_cv"] = out["qc_cv"] > max_cv
    return out


def run_all(ds: Dataset, z_thresh: float = 3.0, ppm: float = 5.0, blank_fold: float = 3.0,
            max_cv: float = 30.0) -> tuple[pd.DataFrame, RTModel]:
    f = annotate(ds)
    rt, model, X = rt_outliers(f, z_thresh)
    parts = [f.drop(columns=["lipid"]), rt, elution_order(f, rt), isobars(f, rt, model, X, ppm, z_thresh=z_thresh),
             duplicates(f, rt, X, model), blank_and_qc(ds, blank_fold, max_cv)]
    res = pd.concat(parts, axis=1)
    flags = [c for c in res if c.startswith("flag_")]
    res[flags] = res[flags].fillna(False).astype(bool)
    res["n_flags"] = res[flags].sum(axis=1)
    res["annotated"] = f["lipid"].notna()
    return res, model
