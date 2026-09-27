"""CSV + self-contained HTML report (figures embedded, opens offline)."""
from __future__ import annotations

import base64
import html
import io
from pathlib import Path

import numpy as np
import pandas as pd

FLAG_TEXT = {
    "flag_rt": "RT off-model",
    "flag_order": "elution order",
    "flag_isobar_rt": "isobar fits RT better",
    "flag_bcfa": "BCFA isomer candidate",
    "flag_blank": "blank contamination",
    "flag_qc_cv": "QC CV high",
    "flag_low_detect": "low detection",
}


def _png(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=130, bbox_inches="tight")
    return base64.b64encode(buf.getvalue()).decode()


def figures(res: pd.DataFrame) -> dict[str, str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = {}
    a = res[res["annotated"] & res["rt_pred"].notna()]
    top = a["lipid_class"].value_counts().index[:10]
    cmap = plt.get_cmap("tab10")

    fig, ax = plt.subplots(figsize=(6, 5.5))
    for i, cls in enumerate(top):
        s = a[a["lipid_class"] == cls]
        ax.scatter(s["rt_pred"], s["rt"], s=8, alpha=0.6, color=cmap(i), label=cls)
    rest = a[~a["lipid_class"].isin(top)]
    ax.scatter(rest["rt_pred"], rest["rt"], s=8, alpha=0.4, color="#999", label="other")
    fl = a[a["flag_rt"]]
    ax.scatter(fl["rt_pred"], fl["rt"], s=40, facecolors="none", edgecolors="k", linewidths=0.8, label="flagged")
    lim = [np.nanmin(a[["rt", "rt_pred"]].to_numpy()), np.nanmax(a[["rt", "rt_pred"]].to_numpy())]
    ax.plot(lim, lim, "k--", lw=0.8)
    r2 = 1 - np.nanvar(a["rt_resid"]) / np.nanvar(a["rt"])
    ax.set(xlabel="predicted RT (min, out-of-fold)", ylabel="observed RT (min)",
           title=f"RT model: R² = {r2:.3f}, median |error| = {a['rt_resid'].abs().median():.2f} min")
    ax.legend(fontsize=7, markerscale=1.5, ncol=2)
    out["rt_scatter"] = _png(fig)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    data = [a.loc[a["lipid_class"] == c, "rt_resid"].dropna() for c in top]
    ax.boxplot(data, tick_labels=list(top), showfliers=True, flierprops={"markersize": 3})
    ax.axhline(0, color="k", lw=0.6)
    ax.set(ylabel="observed - predicted RT (min)", title="Residuals by class")
    ax.tick_params(axis="x", rotation=45)
    out["resid_box"] = _png(fig)
    plt.close(fig)

    if "qc_cv" in res or "sample_blank_ratio" in res:
        fig, axes = plt.subplots(1, 2, figsize=(9, 3.5))
        if "qc_cv" in res:
            axes[0].hist(res["qc_cv"].dropna().clip(upper=100), bins=40, color="#4a6fa5")
            axes[0].axvline(30, color="r", ls="--")
            axes[0].set(xlabel="QC CV %", title="QC precision")
        if "sample_blank_ratio" in res:
            v = np.log10(res["sample_blank_ratio"].replace([np.inf, 0], np.nan).dropna())
            axes[1].hist(v.clip(-1, 4), bins=40, color="#4a6fa5")
            axes[1].axvline(np.log10(3), color="r", ls="--")
            axes[1].set(xlabel="log10 sample / blank", title="Blank contribution")
        out["qc"] = _png(fig)
        plt.close(fig)
    return out


def write(res: pd.DataFrame, out_dir: str | Path, source: str = "", params: dict | None = None) -> Path:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    flags = [c for c in res if c.startswith("flag_")]
    res.to_csv(out_dir / "lipid_qc_all_features.csv")
    flagged = res[res["n_flags"] > 0].copy()
    flagged["reasons"] = flagged[flags].apply(lambda r: ", ".join(FLAG_TEXT.get(c, c) for c in flags if r[c]), axis=1)
    cols = ["name", "adduct", "rt", "rt_pred", "rt_z", "reasons", "order_notes", "best_alt", "best_alt_z", "dup_call",
            "sample_blank_ratio", "qc_cv"]
    cols = [c for c in cols if c in flagged]
    flagged = flagged.sort_values(["n_flags", "rt_z"], key=lambda s: s.abs() if s.name == "rt_z" else s, ascending=False)
    flagged[cols].to_csv(out_dir / "lipid_qc_flagged.csv")

    figs = figures(res)
    counts = res[flags].sum().rename(FLAG_TEXT)
    a = res[res["annotated"]]
    summary = pd.DataFrame({
        "value": [len(res), int(res["annotated"].sum()), a["lipid_class"].nunique(),
                  f"{a['rt_resid'].abs().median():.3f} min", int((res["n_flags"] > 0).sum())]},
        index=["features", "parsed lipid annotations", "lipid classes", "median |RT error|", "features with ≥1 flag"])
    table = flagged[cols].head(300).copy()
    for c in ("rt", "rt_pred", "rt_z", "best_alt_z", "sample_blank_ratio", "qc_cv"):
        if c in table:
            table[c] = pd.to_numeric(table[c], errors="coerce").round(2)
    img = lambda k: f'<img src="data:image/png;base64,{figs[k]}">' if k in figs else ""
    doc = f"""<!doctype html><html><head><meta charset="utf-8"><title>Lipid annotation QC</title>
<style>
body{{font:14px/1.45 -apple-system,Segoe UI,Helvetica,Arial,sans-serif;margin:24px auto;max-width:1200px;padding:0 16px;color:#222;background:#fff}}
h1{{font-size:22px}} h2{{font-size:17px;margin-top:28px}} img{{max-width:100%;margin:6px 8px 6px 0}}
table{{border-collapse:collapse;font-size:12px;width:100%}} th,td{{border-bottom:1px solid #ddd;padding:4px 6px;text-align:left;vertical-align:top}}
th{{background:#f4f6f8;position:sticky;top:0}} .wrap{{overflow-x:auto}} .muted{{color:#666}}
</style></head><body>
<h1>Lipid annotation QC</h1>
<p class="muted">{html.escape(source)} · parameters: {html.escape(str(params or {}))}</p>
<h2>Summary</h2>{summary.to_html(header=False)}
<h2>Flags</h2>{counts.to_frame("features").to_html()}
<h2>Retention-time model</h2>{img("rt_scatter")}{img("resid_box")}
{"<h2>Blank and QC</h2>" + img("qc") if "qc" in figs else ""}
<h2>Flagged features (top 300; full list in lipid_qc_flagged.csv)</h2>
<div class="wrap">{table.to_html(na_rep="", escape=True)}</div>
<h2>How to read this</h2>
<ul>
<li><b>RT off-model</b>: observed RT is &gt;3 robust SD (and &gt;0.2 min) from what this run's own lipids predict for that class/C/DB. Check the peak and MS/MS.</li>
<li><b>Elution order</b>: within a class, more carbons should elute later and more double bonds earlier (reverse phase); Cer should elute before its DHCer. The feature that fits worse is blamed.</li>
<li><b>Isobar fits RT better</b>: another class/adduct with the same m/z (e.g. [M+Na]<sup>+</sup> vs [M+H]<sup>+</sup>, Cer [M+HAc-H]<sup>-</sup> vs hydroxy-Cer [M-H]<sup>-</sup>) matches the observed RT; the current annotation does not.</li>
<li><b>BCFA isomer candidate</b>: the same annotation also appears ~0.3-1 carbon-equivalents earlier, the elution signature of monomethyl branched-chain (iso/anteiso) acyl isomers.</li>
<li>RT cannot separate co-eluting isobars such as PC x:y vs PE (x+3):y; those need MS/MS head-group fragments.</li>
</ul></body></html>"""
    path = out_dir / "lipid_qc_report.html"
    path.write_text(doc)
    return path
