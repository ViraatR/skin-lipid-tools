"""python -m insilico_pace [download|run] [--no-linex]"""
from __future__ import annotations

import argparse
import warnings

import numpy as np
import pandas as pd

from . import data, lipidome, rank, transcriptome
from .catalog import build_catalog


def run(use_linex: bool = True) -> None:
    out = data.RESULTS
    out.mkdir(parents=True, exist_ok=True)

    print("[1/5] enzyme catalog")
    cat = build_catalog()
    cat.drop(columns=["rhea_ids"]).to_csv(out / "enzyme_catalog.csv")

    print("[2/5] transcriptomes (AD/PSO RNA-seq, cytokine KCs, acne, HPA skin)")
    expr = transcriptome.build_expression()

    print("[3/5] lipidome activity proxies (ST002195 Spade vs WT)")
    feats, x = lipidome.load_spade()
    design = lipidome.sample_design(x.index)
    proxies = lipidome.compute_proxies(feats, x)
    effects = lipidome.genotype_effects(proxies, design)
    effects.to_csv(out / "lipid_activity_proxies.csv")
    proxies.join(design).to_csv(out / "lipid_proxy_values_per_animal.csv")
    lipidome.species_effects(feats, x, design).to_csv(out / "lipid_species_spade_vs_wt.csv")
    lip_ev = lipidome.enzyme_lipid_evidence(effects)

    linex = None
    if use_linex:
        print("[4/5] LINEX2 network enrichment")
        r2g: dict[str, set[str]] = {}
        for g, ids in cat["rhea_ids"].items():
            for i in ids:
                r2g.setdefault(i, set()).add(g)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                _, sub, _, pvals = lipidome.linex_enrichment(feats, x, design)
            linex = lipidome.linex_enzymes(sub, pvals, r2g)
            linex.to_csv(out / "linex2_enzymes.csv")
        except Exception as e:  # LINEX2 is optional; the proxy method stands alone
            print(f"  LINEX2 skipped: {e}")

    print("[5/5] integrate + rank")
    ranked = rank.rank_enzymes(cat, expr, lip_ev, linex)
    ranked[rank.REPORT_COLS].to_csv(out / "ranked_candidates.csv")
    sh = ranked[ranked["serine_hydrolase"] & ranked["expressed_in_skin"]]
    sh[rank.REPORT_COLS].to_csv(out / "pace_serine_hydrolase_predictions.csv")
    ctl, p = rank.control_recovery(ranked)
    ctl.to_csv(out / "positive_controls.csv")
    _figures(ranked, effects, proxies, design, out)
    _summary(ranked, sh, ctl, p, effects, linex, expr, out)
    print(f"done -> {out}")


def _figures(ranked, effects, proxies, design, out):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # proxy trajectories: Spade - WT difference by age
    ages = ["P0", "4w", "8w", "10w"]
    top = effects.dropna(subset=["p"]).sort_values("p").head(9).index
    fig, axes = plt.subplots(3, 3, figsize=(10, 8), sharex=True)
    for ax, name in zip(axes.flat, top):
        for gt, col in (("WT", "#4a6fa5"), ("sp", "#c0504d")):
            m = [proxies.loc[design.index[(design.age == a) & (design.genotype == gt)], name] for a in ages]
            ax.errorbar(range(4), [v.mean() for v in m], yerr=[v.std() for v in m], color=col,
                        marker="o", capsize=2, label="Spade" if gt == "sp" else "WT")
        ax.set_title(f"{name}\nq={effects.at[name, 'q']:.1g}", fontsize=9)
        ax.set_xticks(range(4), ages)
    axes.flat[0].legend(fontsize=8)
    fig.supylabel("log2 product / substrate")
    fig.tight_layout()
    fig.savefig(out / "fig_lipid_proxies.png", dpi=150)
    plt.close(fig)

    # top candidates heatmap
    cols = ["AD_LS_vs_HC_lfc", "AD_LS_vs_NL_lfc", "PSO_LS_vs_HC_lfc", "KC_IL4_lfc", "KC_IL13_lfc",
            "KC_IFNg_lfc", "KC_IL17A_lfc", "ACNE_LS_vs_NL_lfc"]
    t = ranked.head(30)
    fig, ax = plt.subplots(figsize=(7, 9))
    im = ax.imshow(t[cols].to_numpy(float), cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
    ax.set_yticks(range(len(t)), [f"{g}{' *' if s else ''}" for g, s in zip(t.index, t["serine_hydrolase"])], fontsize=8)
    ax.set_xticks(range(len(cols)), [c.replace("_lfc", "") for c in cols], rotation=60, ha="right", fontsize=8)
    fig.colorbar(im, ax=ax, label="log2 fold change", shrink=0.5)
    ax.set_title("Top 30 candidates (* = FP-reactive serine hydrolase)", fontsize=10)
    fig.tight_layout()
    fig.savefig(out / "fig_top_candidates.png", dpi=150)
    plt.close(fig)


def _fmt(df: pd.DataFrame) -> str:
    return df.to_markdown(floatfmt=".3g") if hasattr(df, "to_markdown") else df.to_string()


def _summary(ranked, sh, ctl, p, effects, linex, expr, out):
    cols = ["rank", "score", "serine_hydrolase", "ad_direction", "AD_LS_vs_HC_lfc", "th2_lfc",
            "lipid_proxy", "activity_t", "top_skin_cell_type"]
    lines = [
        "# In silico PACE: results summary",
        "",
        f"Enzymes scored: {len(ranked)} ({int(ranked['serine_hydrolase'].sum())} FP-reactive serine hydrolases).",
        "",
        "## Positive controls",
        f"Mean percentile {ctl['percentile'].mean():.0f}; one-sided Mann-Whitney vs all other enzymes p = {p:.2g}.",
        "",
        _fmt(ctl.round(3)),
        "",
        "## Top 25 candidates",
        _fmt(ranked[cols].head(25).round(3)),
        "",
        "## Serine hydrolase predictions",
        "Skin-expressed FP-reactive enzymes, ranked. `th2_lfc` is the mean IL-4/IL-13 transcript change in "
        "primary keratinocytes (GSE36287). Transcript changes are a prior, not a measure of enzyme activity.",
        _fmt(sh[cols + ["KC_nCPM"]].head(20).round(3)),
        "",
        "## Lipidome activity proxies (Spade vs WT, age-adjusted, 4-10 wk)",
        _fmt(effects[["lfc_spade", "t", "q", "lfc_P0", "lfc_4w", "lfc_8w", "lfc_10w", "note"]].sort_values("q").round(3)),
    ]
    if linex is not None and len(linex):
        lines += ["", "## LINEX2 enriched subnetwork enzymes", _fmt(linex.sort_values("linex_edges", ascending=False))]
    (out / "SUMMARY.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser(prog="insilico_pace")
    ap.add_argument("cmd", choices=["download", "run", "all"])
    ap.add_argument("--no-linex", action="store_true")
    a = ap.parse_args()
    if a.cmd in ("download", "all"):
        data.download()
    if a.cmd in ("run", "all"):
        run(use_linex=not a.no_linex)


if __name__ == "__main__":
    main()
