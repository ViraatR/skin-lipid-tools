"""Benchmark: inject realistic annotation errors into real datasets and measure detection.

Error types (each applied to a random 5% of annotated features, one type at a time):
  carbon_pm2   right class, total carbons off by 2 (e.g. PC 34:1 -> PC 36:1)
  db_pm1       right class, one double bond off
  class_isobar swapped to a same-m/z class: PC x:y <-> PE (x+3):y
  na_adduct    [M+H]+ read as [M+Na]+ of a species 2C / 3DB smaller (common in TG, PC)
  wrong_peak   annotation kept, RT replaced by that of a random other feature

Detection score = |robust z| of the RT residual (out-of-fold). Reported: ROC AUC,
recall at the default threshold (|z| > 3 and > 0.2 min), and the false-positive
rate on untouched features (which also contains the dataset's own real errors).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from lipidname import parse

from .checks import annotate, rt_outliers
from .io import Dataset


def _mutate(l, kind: str, rng: np.random.Generator):
    if kind == "carbon_pm2":
        return l.lipid_class, l.total_c + rng.choice([-2, 2]), l.total_db
    if kind == "db_pm1":
        d = l.total_db + rng.choice([-1, 1])
        return (l.lipid_class, l.total_c, d) if d >= 0 else None
    if kind == "class_isobar":
        if l.lipid_class == "PC":
            return "PE", l.total_c + 3, l.total_db
        if l.lipid_class == "PE" and l.total_c > 30:
            return "PC", l.total_c - 3, l.total_db
        return None
    if kind == "na_adduct":
        if l.lipid_class in ("PC", "PE", "TG", "DG", "SM") and l.total_db >= 3:
            return l.lipid_class, l.total_c - 2, l.total_db - 3
        return None
    raise ValueError(kind)


KINDS = ["carbon_pm2", "db_pm1", "class_isobar", "na_adduct", "wrong_peak"]


def inject(ds: Dataset, kind: str, frac: float = 0.05, seed: int = 0) -> tuple[Dataset, pd.Index]:
    rng = np.random.default_rng(seed)
    f = ds.features.copy()
    lip = f["name"].map(parse)
    ok = f.index[lip.notna() & f["rt"].notna()]
    order = rng.permutation(ok)
    n_target = max(5, int(frac * len(ok)))
    hit = []
    for fid in order:
        if len(hit) >= n_target:
            break
        if kind == "wrong_peak":
            other = rng.choice(ok)
            if abs(f.at[other, "rt"] - f.at[fid, "rt"]) < 0.5:
                continue
            f.at[fid, "rt"] = f.at[other, "rt"]
            hit.append(fid)
            continue
        m = _mutate(lip[fid], kind, rng)
        if m is None:
            continue
        cls, c, d = m
        if cls.startswith(("Cer", "SM", "HexCer")) or cls.endswith("_EST"):
            continue  # sphingolipid subclass names need chain detail; keep them out
        name = f"{cls[:-2]} {cls[-1]}-{c}:{d}" if cls.endswith(("-O", "-P")) else f"{cls} {c}:{d}"
        f.at[fid, "name"] = name
        hit.append(fid)
    return Dataset(f, ds.intensity, ds.samples, ds.source), pd.Index(hit)


def benchmark(ds: Dataset, kinds=KINDS, frac: float = 0.05, seeds=(0, 1, 2), z_thresh: float = 3.0) -> pd.DataFrame:
    rows = []
    for kind in kinds:
        for seed in seeds:
            bad_ds, hit = inject(ds, kind, frac, seed)
            f = annotate(bad_ds)
            rt, _, _ = rt_outliers(f, z_thresh)
            scored = rt.dropna(subset=["rt_z"])
            y = scored.index.isin(hit)
            if y.sum() == 0 or y.all():
                continue
            s = scored["rt_z"].abs()
            rows.append({
                "error_type": kind, "seed": seed, "n_features": len(scored), "n_injected": int(y.sum()),
                "auc": roc_auc_score(y, s),
                "recall": scored.loc[y, "flag_rt"].mean(),
                "fpr_clean": scored.loc[~y, "flag_rt"].mean(),
            })
    return pd.DataFrame(rows)


def summarise(bench: pd.DataFrame) -> pd.DataFrame:
    return bench.groupby("error_type")[["n_injected", "auc", "recall", "fpr_clean"]].mean().reindex(
        [k for k in KINDS if k in set(bench["error_type"])])
