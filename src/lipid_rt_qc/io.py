"""Readers: MS-DIAL alignment exports, Skyline reports, and plain CSV.

Everything is normalised to a `Dataset`:
  features  DataFrame indexed by feature_id with name, rt (min), mz, adduct
  intensity DataFrame features x samples (may be empty for annotation-only QC)
  samples   DataFrame indexed by sample with a `type` column: sample | blank | qc | standard
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass
class Dataset:
    features: pd.DataFrame
    intensity: pd.DataFrame = field(default_factory=pd.DataFrame)
    samples: pd.DataFrame = field(default_factory=lambda: pd.DataFrame(columns=["type"]))
    source: str = ""


def _sample_type(name: str, blank_pat: str, qc_pat: str) -> str:
    if re.search(blank_pat, name, re.I):
        return "blank"
    if re.search(qc_pat, name, re.I):
        return "qc"
    return "sample"


def _read_table(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix.lower() in (".txt", ".tsv", ".msp") else ","
    return pd.read_csv(path, sep=sep, header=None, dtype=str, keep_default_na=False, encoding_errors="replace")


# ------------------------------------------------------------------ MS-DIAL
MSDIAL_META = {"Class", "File type", "Injection order", "Batch ID"}
MSDIAL_FILE_TYPES = {"sample": "sample", "blank": "blank", "qc": "qc", "standard": "standard"}


def read_msdial(path: str | Path, blank_pat: str = r"blank|blk|solvent", qc_pat: str = r"\bqc|pool") -> Dataset:
    """MS-DIAL 4/5 alignment result export (Height or Area, .txt)."""
    path = Path(path)
    raw = _read_table(path)
    hdr_row = next((i for i in range(min(len(raw), 20)) if "Alignment ID" in raw.iloc[i].tolist()), None)
    if hdr_row is None:
        raise ValueError(f"{path.name}: no 'Alignment ID' header in the first 20 rows; "
                         "is this an MS-DIAL alignment export?")
    header = raw.iloc[hdr_row].tolist()
    body = raw.iloc[hdr_row + 1:].reset_index(drop=True)
    body.columns = header

    # The metadata block (Class / File type / ...) sits above the header row; its
    # labels are in the column just before the first sample column.
    meta_rows = raw.iloc[:hdr_row]
    label_col = None
    for j in range(raw.shape[1]):
        if set(meta_rows.iloc[:, j]) & MSDIAL_META:
            label_col = j
            break
    sample_cols: list[str] = []
    types: dict[str, str] = {}
    if label_col is not None:
        labels = meta_rows.iloc[:, label_col].tolist()
        ft_row = labels.index("File type") if "File type" in labels else None
        cls_row = labels.index("Class") if "Class" in labels else None
        for j in range(label_col + 1, raw.shape[1]):
            name = header[j]
            ft = meta_rows.iat[ft_row, j].strip().lower() if ft_row is not None else ""
            cls = meta_rows.iat[cls_row, j].strip() if cls_row is not None else ""
            if not name or (not ft and not cls):
                break  # past the samples, into Average/Stdev summary columns
            sample_cols.append(name)
            types[name] = MSDIAL_FILE_TYPES.get(ft, _sample_type(name, blank_pat, qc_pat))

    rt_col = next((c for c in header if c.startswith("Average Rt")), None)
    mz_col = next((c for c in header if c.startswith("Average Mz")), None)
    if rt_col is None or mz_col is None:
        raise ValueError(f"{path.name}: MS-DIAL export needs 'Average Rt(min)' and 'Average Mz' columns")
    feats = pd.DataFrame({
        "feature_id": body["Alignment ID"],
        "name": body["Metabolite name"],
        "rt": pd.to_numeric(body[rt_col], errors="coerce"),
        "mz": pd.to_numeric(body[mz_col], errors="coerce"),
        "adduct": body.get("Adduct type", ""),
        "ontology": body.get("Ontology", ""),
        "msms_matched": body.get("MS/MS matched", ""),
    }).set_index("feature_id")
    inten = body[sample_cols].apply(pd.to_numeric, errors="coerce")
    inten.index = feats.index
    samples = pd.DataFrame({"type": pd.Series(types)}).reindex(sample_cols)
    return Dataset(feats, inten, samples, source=f"MS-DIAL: {path.name}")


# ------------------------------------------------------------------ Skyline
SKY_NAME = ["Molecule Name", "Molecule", "Peptide", "Precursor Name"]
SKY_RT = ["Retention Time", "Best Retention Time", "Peak Rt", "Explicit Retention Time"]
SKY_MZ = ["Precursor Mz", "Precursor m/z", "Mz"]
SKY_AREA = ["Area", "Total Area", "Total Area MS1", "Peak Area"]
SKY_REP = ["Replicate Name", "Replicate", "File Name"]
SKY_ADDUCT = ["Precursor Adduct", "Adduct"]


def _pick(cols, options):
    return next((c for c in options if c in cols), None)


def read_skyline(path: str | Path, blank_pat: str = r"blank|blk|solvent", qc_pat: str = r"\bqc|pool") -> Dataset:
    """Skyline custom report (long format: one row per molecule x replicate)."""
    path = Path(path)
    df = pd.read_csv(path, sep=None, engine="python", encoding_errors="replace")
    cols = list(df.columns)
    name, rt, mzc = _pick(cols, SKY_NAME), _pick(cols, SKY_RT), _pick(cols, SKY_MZ)
    area, rep, add = _pick(cols, SKY_AREA), _pick(cols, SKY_REP), _pick(cols, SKY_ADDUCT)
    if name is None or rt is None:
        raise ValueError(f"Skyline report needs a molecule name and retention time column; got {cols}")
    df[rt] = pd.to_numeric(df[rt], errors="coerce")
    if area:
        df[area] = pd.to_numeric(df[area], errors="coerce")
    key = [name] + ([add] if add else [])
    g = df.groupby(key, dropna=False)
    feats = g.agg(rt=(rt, "median"), **({"mz": (mzc, "first")} if mzc else {})).reset_index()
    feats = feats.rename(columns={name: "name", **({add: "adduct"} if add else {})})
    if "adduct" not in feats:
        feats["adduct"] = ""
    if "mz" not in feats:
        feats["mz"] = np.nan
    feats["feature_id"] = [f"sky{i:05d}" for i in range(len(feats))]
    feats = feats.set_index("feature_id")
    inten, samples = pd.DataFrame(index=feats.index), pd.DataFrame(columns=["type"])
    if area and rep:
        wide = df.pivot_table(index=key, columns=rep, values=area, aggfunc="sum")
        wide.index = wide.index.map(lambda k: k if isinstance(k, tuple) else (k,))
        lookup = {tuple(r[["name"] + (["adduct"] if add else [])]): fid for fid, r in feats.iterrows()}
        wide.index = [lookup.get(k) for k in wide.index]
        inten = wide.reindex(feats.index)
        samples = pd.DataFrame({"type": [_sample_type(c, blank_pat, qc_pat) for c in inten.columns]}, index=inten.columns)
    return Dataset(feats, inten, samples, source=f"Skyline: {path.name}")


# ------------------------------------------------------------------ generic
def read_generic(path: str | Path, name_col="name", rt_col="rt", mz_col="mz", adduct_col="adduct",
                 blank_pat: str = r"blank|blk|solvent", qc_pat: str = r"\bqc|pool") -> Dataset:
    """Wide CSV: annotation columns + one numeric column per sample."""
    path = Path(path)
    df = pd.read_csv(path, sep=None, engine="python")
    missing = [c for c in (name_col, rt_col) if c not in df]
    if missing:
        raise ValueError(f"{path.name}: CSV input needs {missing} column(s); got {list(df.columns)}")
    feats = pd.DataFrame({
        "name": df[name_col], "rt": pd.to_numeric(df[rt_col], errors="coerce"),
        "mz": pd.to_numeric(df[mz_col], errors="coerce") if mz_col in df else np.nan,
        "adduct": df[adduct_col] if adduct_col in df else "",
    })
    feats.index = [f"f{i:05d}" for i in range(len(df))]
    feats.index.name = "feature_id"
    meta = {name_col, rt_col, mz_col, adduct_col}
    num = [c for c in df.columns if c not in meta and pd.api.types.is_numeric_dtype(df[c])]
    inten = df[num].set_axis(feats.index)
    samples = pd.DataFrame({"type": [_sample_type(c, blank_pat, qc_pat) for c in num]}, index=num)
    return Dataset(feats, inten, samples, source=f"CSV: {path.name}")


# ------------------------------------------------------------------ MetaboLights MAF
def read_maf(path: str | Path, blank_pat: str = r"blank|blk|solvent", qc_pat: str = r"\bqc|pool") -> Dataset:
    """MetaboLights metabolite assignment file (m_*.tsv). Sample columns follow the
    `smallmolecule_abundance_*` block; the adduct is usually in `modifications`."""
    path = Path(path)
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = [c for c in ("metabolite_identification", "retention_time", "mass_to_charge") if c not in df]
    if missing:
        raise ValueError(f"{path.name}: MetaboLights MAF needs {missing} column(s); got {list(df.columns)}")
    adduct = pd.Series("", index=df.index)
    for col in ("modifications", "fragmentation", "adduct"):
        if col in df:
            v = df[col].where(df[col].str.match(r"^\[M.*\][+-]$"), "")
            adduct = adduct.where(adduct != "", v)
    feats = pd.DataFrame({
        "name": df["metabolite_identification"],
        "rt": pd.to_numeric(df["retention_time"], errors="coerce"),
        "mz": pd.to_numeric(df["mass_to_charge"], errors="coerce"),
        "adduct": adduct,
    })
    if feats["rt"].max() > 120:  # submitted in seconds
        feats["rt"] = feats["rt"] / 60
    feats.index = pd.Index([f"maf{i:05d}" for i in range(len(df))], name="feature_id")
    abund = [i for i, c in enumerate(df.columns) if c.startswith("smallmolecule_abundance")]
    num = list(df.columns[max(abund) + 1:]) if abund else []
    inten = df[num].apply(pd.to_numeric, errors="coerce").set_axis(feats.index)
    samples = pd.DataFrame({"type": [_sample_type(c, blank_pat, qc_pat) for c in num]}, index=num)
    return Dataset(feats, inten, samples, source=f"MetaboLights MAF: {path.name}")


def read_any(path: str | Path, fmt: str = "auto", **kw) -> Dataset:
    path = Path(path)
    if fmt == "auto":
        head = path.read_text(errors="replace")[:20000]
        fmt = ("msdial" if "Alignment ID" in head else "maf" if "metabolite_identification" in head
               else "skyline" if any(c in head for c in SKY_NAME[:2] + SKY_REP) else "csv")
    return {"msdial": read_msdial, "maf": read_maf, "skyline": read_skyline, "csv": read_generic}[fmt](path, **kw)
